#!/usr/bin/env python
from src.group import readBam, read_bam_from_bed
from src.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from src.get_consensus import get_cons_dict, get_all_consensus, write_singleton_reads, get_reference_sequence
from src.get_cons_info import get_cons_info, write_consensus, calc_major_nonref_allele_frequency
from src.get_regions_from_bed import read_bed, sort_regions, merge_regions, get_overlap
import sys
import os
import pysam
from multiprocessing import Pool, cpu_count
import argparse
import logging

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)")
    parser.add_argument('-o',  '--output_path', dest='output_path', 
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-b', '--bam', dest='bam_file', help='Path to BAM-file')
    parser.add_argument('-bed', '--bed-file', dest='bed_file', 
                        help='Path to a BED file defining the targeted regions, \
                              i.e. chromosomal positions. The Bed file is used \
                              for annotation.')
    parser.add_argument('-regions_from_bed', dest='regions_from_bed', 
                        help='Include this flag if regions used for UMI clustering \
                              and variant calling should be defined from the BED file. \
                              Default is to detect the regions automatically from the \
                              BAM file. ', 
                        action='store_true')
    parser.add_argument('-r', '--reference', dest='reference_file', 
                        help='Path to the reference sequence in Fasta format, Used for annotation, required', required=True)
    parser.add_argument('-s', '--sample_name', dest='sample_name', 
                        help='Sample name that will be used as base name for the output files. \
                        If excluded the sample name will be extracted from the BAM file.')
    parser.add_argument('-d', '--edit_distance', dest='edit_distance_threshold', 
                        help="Edit distance threshold for UMI clustering, [default = %(default)s]", default=1)
    parser.add_argument('-p', '--position_threshold', dest='position_threshold', 
                        help='Position threshold for grouping by position [default = %(default)s]', default=10)
    parser.add_argument('-cons_freq', '--consensus_frequency_threshold', dest='consensus_frequency_threshold', 
                        help='Minimum percent of the majority base at a position for consensus to be called. \
                              [default = %(default)s]', default=60.0)

    parser.add_argument('-indel_freq', '--indel_frequency_threshold', dest='indel_frequency_threshold', 
                        help='Percent threshold for indels to be included in the consensus read. \
                        [default = %(default)s]', default=60.0)
    parser.add_argument('-singletons', '--include_singletons', dest='include_singletons', action='store_true', 
                        help='Include this flag if singleton reads should be included in the output consensus \
                              read bam file. Note that the singletons will not be error corrected')
    parser.add_argument('-t', '--num_threads', dest='num_threads', 
                        help='Number of threads to run the program on. If excluded, the number of cpus are \
                        automatically detected')
    args = parser.parse_args(sys.argv[1:])

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMI Error Correct')

    return(args)


def check_output_directory(outdir):
    '''Check if outdir exists, otherwise create it'''
    if os.path.isdir(outdir):
        return(outdir)
    else:
        os.mkdir(outdir)
        return(outdir)


def get_sample_name(bamfile):
    '''Get the sample name as the basename of the input files.'''
    sample_name=bamfile.split('/')[-1]    
    if '.sorted' in sample_name:
        sample_name = sample_name.replace('.sorted','')
    sample_name = sample_name.replace('.bam','')
    return(sample_name)


def cluster_consensus_worker(args):
    '''Run UMI clustering and consensus read generation on one region'''
    umi_dict, samplename, tmpfilename, regionid, contig, start, end, edit_distance_threshold, \
    bamfilename, include_singletons, annotations, fasta, indel_frequency_cutoff, \
    consensus_frequency_cutoff = args  # extract args
    
    #UMI clustering
    adj_matrix = cluster_barcodes(umi_dict, edit_distance_threshold)
    clusters = get_connected_components(umi_dict, adj_matrix)
    umis = merge_clusters(umi_dict, clusters)
    
    #Consensus sequence generation
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, contig, 
                                                      start, end, True)  # include_singletons=True
    consensus_seq = get_all_consensus(position_matrix, umis, 
                                      contig, regionid, indel_frequency_cutoff, 
                                      consensus_frequency_cutoff)
    outfilename = tmpfilename
    
    #Write consensus reads (and singletons) to a BAM file
    num_cons=0
    with pysam.AlignmentFile(bamfilename, 'rb') as f, pysam.AlignmentFile(outfilename, 'wb', template=f) as g:
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
                num_cons+=1
        if include_singletons:
            write_singleton_reads(singleton_matrix, regionid, g)
    
    #Generate info for cons file
    cons = get_cons_info(consensus_seq, singleton_matrix)
    consfilename = outfilename.rstrip('.bam') + '.cons'
    statfilename = outfilename.rstrip('.bam') + '.hist'
    if len(cons) > 0:
        startpos = min(list(cons.keys())) #take the rightmost coordinate as start
        endpos = max(list(cons.keys())) + 1 #take the leftmost coordinate as end
        with pysam.FastaFile(fasta) as f:
            ref_seq = get_reference_sequence(f, contig, startpos, endpos)
        with open(consfilename, 'w') as g:
            write_consensus(g, cons, ref_seq, startpos, contig, annotations, samplename, False)
    else:  # empty file
        g = open(consfilename, 'w')
        g.close()
    #Write to hist/stat file
    if len(cons) > 0:
        with open(statfilename, 'w') as g2:
            name=get_overlap(annotations, start, endpos)
            regionname = '{}:{}-{}'.format(contig, start, endpos)
            g2.write('\t'.join([str(regionid), regionname, name, 'consensus_reads: ' + str(num_cons),  
                                'singletons: ' + str(len(singleton_matrix))]) + '\n')
    else:  # empty file
        g = open(statfilename, 'w')
        g.close()

def update_bam_header(bamfile, samplename):
    with pysam.AlignmentFile(bamfile,'rb') as f:
        new_header = f.header.copy().to_dict()
        template = { 'ID': 'L1', 
                     'SM': samplename,
                     'LB': samplename,
                     'PL': 'ILLUMINA'}
    
    new_header['RG'] = [template]
    return(new_header)


def merge_bams(output_path, bamfilelist, sample_name):
    '''Merge all BAM files for in bamfilelist, and remove temporary files'''
    new_header = update_bam_header(bamfilelist[0], sample_name)
    with pysam.AlignmentFile(output_path + '/' + sample_name + '_consensus_reads.bam', 'wb', header=new_header) as g:
        for filename in bamfilelist:
            with pysam.AlignmentFile(filename, 'rb') as f1:
                for line in f1:
                    g.write(line)

    for filename in bamfilelist:
        os.remove(filename)

def merge_cons(output_path, consfilelist, sample_name):
    '''Merge all cons files in consfilelist and remove temporary files.'''
    with open(output_path + '/' + sample_name + '.cons', 'w') as g:
        g.write('Sample Name\tContig\tPosition\tName\tReference\tA\tC\tG\tT\tI\tD\tN\tCoverage\tConsensus group size\tMax Non-ref Allele Count\tMax Non-ref Allele Frequency\tMax Non-ref Allele\n')
        for filename in consfilelist:
            with open(filename) as f:
                for line in f:
                    g.write(line)

    for filename in consfilelist:
        os.remove(filename)

def check_duplicate_positions(cons_file):
    with open(cons_file) as f:
        line = f.readline()
        a = 13
        b = 2
        positions = []
        duppos = []
        for line in f:
            parts = line.split('\t')
            if parts[a] == '0':
                if parts[b] in positions:
                    duppos.append(parts[b])
                positions.append(parts[b])    
    return(duppos)

def sum_lists(*args):
    return(list(map(sum, zip(*args))))

def merge_duplicate_positions(duppos,cons_file):
    dupcons={}
    a = 13
    b = 2
    with open(cons_file) as f:
        line = f.readline()
        for line in f:
            parts = line.split('\t')
            pos = parts[b]
            if pos in duppos:
                fsize = parts[a]
                if pos not in dupcons:
                    dupcons[pos] = {}
                if fsize not in dupcons[pos]:
                    dupcons[pos][fsize]=[]
                dupcons[pos][fsize].append(line)
    #print(dupcons)
    newpos={}
    for pos in dupcons:
        newpos[pos]={}
        for fsize in dupcons[pos]:
            newpos[pos][fsize]=[0,0,0,0,0,0,0,0,0,0]
            for s in dupcons[pos][fsize]:
                parts=s.split('\t')
                newpos[pos][fsize]=[int(a)+int(b) for a,b in zip(newpos[pos][fsize],parts[5:15])]
    with open(cons_file) as f,open(cons_file+'_new', 'w') as g:
        line = f.readline()
        g.write(line)
        positions=[]
        fsizes=['0', '1', '2', '3', '4', '5', '7', '10', '20', '30']
        for line in f:
            parts=line.split('\t')
            pos = parts[b]
            if pos not in newpos:
                g.write(line)
            else:
                if pos not in positions:
                    for fsize in fsizes:
                        if fsize in newpos[pos]:
                            tmp = newpos[pos][fsize]
                            newlist = [str(x) for x in tmp]
                            consdict = { 'A': tmp[0], 'C': tmp[1], 'G':tmp[2], 'T': tmp[3], 'I': tmp[4], 'D': tmp[5], 'N': tmp[6]}
                            tot = sum(consdict.values())
                            mna, freq, count = calc_major_nonref_allele_frequency(consdict, parts[4], tot)
                        #frac=(newpos[pos][fsize][9]/newpos[pos][fsize][7])*1.0
                            g.write('\t'.join(parts[0:5])+'\t'+'\t'.join(newlist[0:8])+'\t'+fsize+'\t'+str(count)+'\t'+str(freq)+'\t'+mna+'\n')
                    positions.append(pos)
    
    os.remove(cons_file)
    os.rename(cons_file+'_new',cons_file)

def merge_stat(output_path, statfilelist, sample_name):
    '''Merge all stat files in statfilelist and remove temporary files.'''
    with open(output_path + '/' + sample_name + '.hist', 'w') as g:
        for filename in statfilelist:
            with open(filename) as f:
                for line in f:
                    g.write(line)

    for filename in statfilelist:
        os.remove(filename)


def index_bam_file(filename, num_threads=1):
    '''Index the consensus reads bam file'''
    pysam.sort('-@',  num_threads, filename, '-o', filename + '.sorted', catch_stdout=False)
    os.rename(filename + '.sorted', filename)
    pysam.index(filename, catch_stdout=False)


def cluster_umis_all_regions(regions, ends, edit_distance_threshold, samplename,  bamfilename, output_path, 
                             include_singletons, fasta, bedregions, num_cpus, 
                             indel_frequency_cutoff, consensus_frequency_cutoff):
    '''Function for running UMI cluestering and error correction using num_cpus threads,
        i.e. one region on each thread.'''
    argvec = []
    bamfilelist = []
    i = 0
    for contig in regions:
        for pos in regions[contig]:
            if contig in bedregions:
                annotations = bedregions[contig]
            else:
                annotations = []
            
            tmpfilename = '{}/tmp_{}.bam'.format(output_path, i)
            argvec.append((regions[contig][pos], samplename, tmpfilename, int(i), contig, int(pos), 
                           int(ends[contig][pos]), int(edit_distance_threshold), bamfilename,
                           include_singletons, annotations, fasta, indel_frequency_cutoff,
                           consensus_frequency_cutoff))
            bamfilelist.append('{}/tmp_{}.bam'.format(output_path, i))
            i += 1

    p = Pool(int(num_cpus))
    
    p.map(cluster_consensus_worker, argvec)
    return(bamfilelist)


def cluster_umis_on_position(bamfilename, position_threshold, group_method, bedfilename=None):
    '''Function for cluster umis on position'''
    position_threshold = 20
    # group_method='fromBed'
    # group_method='automatic'
    if group_method == 'fromBed':
        regions, ends = read_bam_from_bed(bamfilename, bedfilename, position_threshold)
    else:
        regions, ends = readBam(bamfilename, position_threshold)
    
    return(regions, ends)


def run_umi_errorcorrect(args):
    '''Run the umi clustering and consensus read generation (error correction)'''
    logging.info("Starting UMI clustering")    
    args.output_path = check_output_directory(args.output_path)
    if args.regions_from_bed:
        group_method = 'fromBed'
    else:
        group_method = 'automatic'

    logging.info('Group by position method: {}'.format(group_method))
    if not args.sample_name:
        args.sample_name = get_sample_name(args.bam_file)
    regions, ends = cluster_umis_on_position(args.bam_file, args.position_threshold, 
                                             group_method, args.bed_file)
    nregions = 0
    for chrx in regions:
        nregions += len(regions[chrx])
    logging.info("Number of regions, {}".format(nregions))
    
    edit_distance_threshold = args.edit_distance_threshold
    if args.num_threads:
        num_cpus = int(args.num_threads)
    else:
        num_cpus = int(cpu_count())
    logging.info("Starting Consensus sequence generation")
    logging.info("Starting {} threads".format(num_cpus))
    fasta = args.reference_file
    if args.bed_file:
        bedregions = read_bed(args.bed_file)
        bedregions = sort_regions(bedregions)
        if group_method=='fromBed':
            bedregions = merge_regions(bedregions, 0)
    else:
        bedregions = []
    bamfilelist = cluster_umis_all_regions(regions, ends, edit_distance_threshold, 
                                           args.sample_name, args.bam_file, args.output_path, 
                                           args.include_singletons, fasta, bedregions, 
                                           num_cpus, args.indel_frequency_threshold, 
                                           args.consensus_frequency_threshold)
    merge_bams(args.output_path, bamfilelist, args.sample_name)
    index_bam_file(args.output_path + '/' + args.sample_name + '_consensus_reads.bam',
              args.num_threads)
    consfilelist = [x.rstrip('.bam') + '.cons' for x in bamfilelist]
    merge_cons(args.output_path, consfilelist, args.sample_name)
    cons_file = args.output_path + '/' + args.sample_name + '.cons'
    #duppos = check_duplicate_positions(cons_file)
    #if any(duppos):
    #    merge_duplicate_positions(duppos,cons_file)

    statfilelist = [x.rstrip('.bam') + '.hist' for x in bamfilelist]
    merge_stat(args.output_path, statfilelist, args.sample_name)
    logging.info("Consensus generation complete, output written to {}, {}".format(args.output_path + 
                 '/' + args.sample_name + '_consensus_reads.bam',
                 args.output_path + '/' + args.sample_name + '.cons'))

def main(args):
    run_umi_errorcorrect(args)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
