#!/usr/bin/env python
from umierrorcorrect.src.group import readBam, read_bam_from_bed, read_bam_from_tag
from umierrorcorrect.src.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from umierrorcorrect.src.get_consensus3 import get_cons_dict, get_all_consensus, write_singleton_reads, get_reference_sequence
from umierrorcorrect.src.get_cons_info import get_cons_info, write_consensus, calc_major_nonref_allele_frequency
from umierrorcorrect.src.get_regions_from_bed import read_bed, sort_regions, merge_regions, get_overlap
from umierrorcorrect.version import __version__
import sys
import os
import re
import pysam
from multiprocessing import Pool, cpu_count
import subprocess
import argparse
import logging
import pickle
from itertools import islice

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="UmiErrorCorrect v. {}. \
                                                  Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)".format(__version__))
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
    parser.add_argument('-regions_from_tag', dest='regions_from_tag',
                        help='Include this flag if regions used for UMI clustering \
                              and variant calling should be defined from the UG tag in the BAM file. \
                              Default is to detect the regions automatically from the \
                              BAM file. ',
                        action='store_true')
    parser.add_argument('-r', '--reference', dest='reference_file', 
                        help='Path to the reference sequence in Fasta format, Used for annotation, required', required=True)
    parser.add_argument('-s', '--sample_name', dest='sample_name', 
                        help='Sample name that will be used as base name for the output files. \
                        If excluded the sample name will be extracted from the BAM file.')
    parser.add_argument('-remove', '--remove_large_files',  dest='remove_large_files', action='store_true',\
                        help='Include this flag to emove the original Fastq and BAM files (reads without error correction).')
    parser.add_argument('-d', '--edit_distance', dest='edit_distance_threshold', 
                        help="Edit distance threshold for UMI clustering, [default = %(default)s]", default=1)
    parser.add_argument('-p', '--position_threshold', dest='position_threshold', 
                        help='Position threshold for grouping by position [default = %(default)s]', default=20)
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
    
    indel_frequency_cutoff = float(indel_frequency_cutoff)
    consensus_frequency_cutoff = float(consensus_frequency_cutoff)
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
    with open(output_path + '/' + sample_name + '_cons.tsv', 'w') as g:
        g.write('Sample Name\tContig\tPosition\tName\tReference\tA\tC\tG\tT\tI\tD\tN\tCoverage\tConsensus group size\tMax Non-ref Allele Count\tMax Non-ref Allele Frequency\tMax Non-ref Allele\n')
        for filename in consfilelist:
            with open(filename) as f:
                for line in f:
                    g.write(line)

    for filename in consfilelist:
        os.remove(filename)

def check_duplicate_positions(cons_file):
    
    chrlist = []
    with open(cons_file) as f, open('tmp.txt','w') as g:
        f.readline()
        for line in f:
            parts=line.split('\t')
            if parts[1] not in chrlist:
                chrlist.append(parts[1])
            if parts[13]=='0':
                g.write(' '.join(parts[1:3])+'\n')
    command1=['sort tmp.txt | uniq -d']
    with open('tmp2.txt','w') as g:
        p1=subprocess.Popen(command1, shell=True, stdout=g)
        p1.communicate()
    os.remove('tmp.txt')
    duppos={}
    for chrx in chrlist:
        duppos[chrx] = []
    with open('tmp2.txt') as f:
        for line in f:
            line=line.rstrip()
            parts=line.split()
            chrx=parts[0]
            pos=parts[1]
            duppos[chrx].append(pos)            
    os.remove('tmp2.txt')
    return(duppos)

def sum_lists(*args):
    return(list(map(sum, zip(*args))))

def merge_duplicate_positions(args):
    chrx, duppos, cons_file = args
    dupcons={}
    a = 13
    b = 2
    with open(cons_file) as f:
        line = f.readline()
        for line in f:
            parts = line.split('\t')
            pos = parts[b]
            contig = parts[b-1]
            if contig == chrx and pos in duppos:
                fsize = parts[a]
                if pos not in dupcons:
                    dupcons[pos] = {}
                if fsize not in dupcons[pos]:
                    dupcons[pos][fsize]=[]
                dupcons[pos][fsize].append(line)
    newpos={}
    for pos in dupcons:
        newpos[pos]={}
        for fsize in dupcons[pos]:
            newpos[pos][fsize]=[0,0,0,0,0,0,0,0,0,0]
            for s in dupcons[pos][fsize]:
                parts=s.split('\t')
                newpos[pos][fsize]=[int(a)+int(b) for a,b in zip(newpos[pos][fsize],parts[5:15])]
    with open(cons_file) as f,open(cons_file+'_new'+chrx, 'w') as g:
        line = f.readline()
        g.write(line)
        positions=[]
        fsizes=['0', '1', '2', '3', '4', '5', '7', '10', '20', '30']
        for line in f:
            parts=line.split('\t')
            pos = parts[b]
            contig = parts[b-1]
            if contig == chrx:
                if pos not in newpos:
                    g.write(line)
                else:
                    if pos not in positions:
                        for fsize in fsizes:
                            if fsize in newpos[pos]:
                                tmp = newpos[pos][fsize]
                                newlist = [str(x) for x in tmp]
                                consdict = { 'A': tmp[0], 'C': tmp[1], 'G':tmp[2], 'T': tmp[3], 'I': tmp[4], 'D': tmp[5], 'N': tmp[6]}
                                #tot = sum(consdict.values())
                                tot = sum(consdict[key] for key in consdict if key != 'I')
                                if tot>0:
                                    refbase = parts[4]
                                    nonrefcons = {key: consdict[key] for key in consdict if key != refbase}
                                    mna, freq, count = calc_major_nonref_allele_frequency(nonrefcons, parts[4], tot)
                        #frac=(newpos[pos][fsize][9]/newpos[pos][fsize][7])*1.0
                                    g.write('\t'.join(parts[0:5])+'\t'+'\t'.join(newlist[0:8])+'\t'+fsize+'\t'+str(count)+'\t'+str(freq)+'\t'+mna+'\n')
                        positions.append(pos)

    #os.remove(cons_file)
    #os.rename(cons_file+'_new',cons_file)

def merge_duplicate_positions_all_chromosomes(duppos, cons_file, num_cpus):
    argvec=[]
    print(duppos)
    for chrx in duppos:
        tmpargs=(chrx, duppos[chrx], cons_file)
        argvec.append(tmpargs)
    p = Pool(int(num_cpus))
    p.map(merge_duplicate_positions, argvec)
    merge_tmp_cons_files(duppos.keys(), cons_file)
    if os.path.isfile(cons_file+'2'):
        os.remove(cons_file)
        os.rename(cons_file+'2',cons_file)
        

def merge_tmp_cons_files(chrlist, cons_file):
    try:
        chrlist_sorted=sorted(chrlist, key=int)
    except ValueError as e:
        chrlist_sorted=sorted(chrlist)
    tmpfilelist=[cons_file + '_new' + str(x) for x in chrlist_sorted]
    with open(cons_file + '2', 'w') as g:
        i = 0
        for filename in tmpfilelist:
            with open(filename) as f:
                if i > 0:
                    f.readline()
                for line in f:
                    g.write(line)
                i += 1
    for filename in tmpfilelist:
        os.remove(filename)

def merge_stat(output_path,statfilelist, sample_name):
    '''Merge all stat files in statfilelist and remove temporary files.'''
    with open(output_path + '/' + sample_name + '.hist', 'w') as g:
        for filename in statfilelist:
            with open(filename) as f:
                for line in f:
                    g.write(line)

    for filename in statfilelist:
        os.remove(filename)

def merge_duplicate_stat(output_path,samplename):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    histfile=output_path + '/' + samplename + '.hist'
    regions={}
    #histfile
    with open(histfile) as f:
        for line in f:
            line=line.rstrip()
            parts=line.split('\t')
            idx=parts[0]
            chrx=parts[1].split(':')[0]
            if chrx not in regions:
                regions[chrx]={}
            pos=parts[1].split(':')[1].split('-')[0]
            end=parts[1].split(':')[1].split('-')[1]
            name=parts[2]
            numcons=parts[3].split()[1]
            numsing=parts[4].split()[1]
            if pos not in regions[chrx]:
                regions[chrx][pos]=(idx,int(pos),int(end),name,int(numcons),int(numsing))
            else:
                tmp=regions[chrx][pos]
                if '-' in tmp[0]:
                    newid=tmp[0].split('-')[0]+'-'+idx
                else:
                    newid=tmp[0]+'-'+idx
                if int(end) > int(tmp[2]):
                    newend=end
                else:
                    newend=tmp[2]
                newnumcons=tmp[4] + int(numcons)
                newnumsing=tmp[5] + int(numsing)
                regions[chrx][pos]=(newid, tmp[1], newend, name, newnumcons, newnumsing)
    with open(histfile+'2','w') as g:
        for chrx in sorted(regions,key=alphanum_key):
            for pos in regions[chrx]:
                tmp=regions[chrx][pos]
                g.write('{}\t{}:{}-{}\t{}\tconsensus_reads: {}\tsingletons: {}\n'.format(tmp[0],str(chrx),tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]))
    os.rename(histfile+'2', histfile)


def index_bam_file(filename, num_threads=1):
    '''Index the consensus reads bam file'''
    pysam.sort('-@', str(num_threads), filename, '-o', filename + '.sorted', catch_stdout=False)
    os.rename(filename + '.sorted', filename)
    pysam.index(filename, catch_stdout=False)

def split_into_chunks(umi_dict,clusters):
    ''' If one region contains more than 100000 raw reads, split in chunks of 100000.
        keep all barcodes in same cluster in the same chunk.
    '''
    n = 0
    i = 0
    newdicts = []
    newdict = {}
    for c in clusters:
        for j in c:
            count=umi_dict[j]
            newdict[j]=count
            n+=count
        if n > 100000:
            newdicts.append(newdict)
            newdict={}
            n=0
    newdicts.append(newdict) #add remaining entries
    return(newdicts)
    
    
    #n = 0
    #i = 0
    #newdicts = []
    #b=list(dictname.values())
    #a=iter(dictname.items())
    #for count in b:
    #    n += count
    #    i += 1
    #    if n > 100000:
    #        newdicts.append(dict(islice(a,i))) #add (more than) 100000 raw reads to newdicts (from 0 to index i)
    #        n = 0
    #        i = 0
    #newdicts.append(dict(islice(a,i))) #add remaining entries
    #return(newdicts)

def cluster_umis_all_regions(regions, ends, edit_distance_threshold, samplename,  bamfilename, output_path, 
                             include_singletons, fasta, bedregions, num_cpus, 
                             indel_frequency_cutoff, consensus_frequency_cutoff, region_from_tag=False,starts=[]):
    '''Function for running UMI cluestering and error correction using num_cpus threads,
        i.e. one region on each thread.'''
    argvec = []
    bamfilelist = []
    i = 0
    j = 0
    for contig in regions:
        for pos in regions[contig]:
            if contig in bedregions:
                annotations = bedregions[contig]
            else:
                annotations = []
            
            if region_from_tag:
                i=pos
                posx=int(starts[contig][pos])
                j=0
            else:
                posx=int(pos)
            tmpfilename = '{}/tmp_{}.bam'.format(output_path, i)
            numreads = sum(regions[contig][pos].values())
            if numreads > 100000: #split in chunks
                umi_dict=regions[contig][pos]
                adj_matrix = cluster_barcodes(umi_dict, edit_distance_threshold)
                clusters = get_connected_components(umi_dict, adj_matrix)
                newdicts=split_into_chunks(umi_dict,clusters)
                for x in newdicts:
                    tmpfilename = '{}/tmp_{}.bam'.format(output_path, i)
                    argvec.append((x, samplename, tmpfilename, i, contig, posx,
                                    int(ends[contig][pos]), int(edit_distance_threshold), bamfilename,
                                    include_singletons, annotations, fasta, indel_frequency_cutoff,
                                    consensus_frequency_cutoff))
                    bamfilelist.append('{}/tmp_{}.bam'.format(output_path, i))
                    if not region_from_tag:
                        i += 1
                    else:
                        i=i+'_'+str(j)
                        j += 1
            else:
                argvec.append((regions[contig][pos], samplename, tmpfilename, i, contig, posx, 
                                int(ends[contig][pos]), int(edit_distance_threshold), bamfilename,
                                include_singletons, annotations, fasta, indel_frequency_cutoff,
                                consensus_frequency_cutoff))
                bamfilelist.append('{}/tmp_{}.bam'.format(output_path, i))
                if not region_from_tag:
                    i += 1

    p = Pool(int(num_cpus))
    
    p.map(cluster_consensus_worker, argvec)
    return(bamfilelist)


def cluster_umis_on_position(bamfilename, position_threshold, group_method, bedfilename=None):
    '''Function for 0cluster umis on position'''
    #position_threshold = 20
    # group_method='fromBed'
    # group_method='automatic'
    position_threshold=int(position_threshold)
    if group_method == 'fromBed':
        regions, ends = read_bam_from_bed(bamfilename, bedfilename, position_threshold)
    elif group_method == 'fromTag':
        regions, starts, ends = read_bam_from_tag(bamfilename)
    else:
        regions, ends = readBam(bamfilename, position_threshold)
    
    if group_method == 'fromTag':
        return(regions,ends,starts)
    else:
        return(regions, ends)


def run_umi_errorcorrect(args):
    '''Run the umi clustering and consensus read generation (error correction)'''
    logging.info("Starting UMI clustering")    
    args.output_path = check_output_directory(args.output_path)
    if args.regions_from_bed:
        group_method = 'fromBed'
    elif args.regions_from_tag:
        group_method = 'fromTag'
    else:
        group_method = 'automatic'

    logging.info('Group by position method: {}'.format(group_method))
    if not args.sample_name:
        args.sample_name = get_sample_name(args.bam_file)
    if group_method == 'fromTag':
        regions, ends, starts = cluster_umis_on_position(args.bam_file, args.position_threshold,
                                             group_method, args.bed_file)
    else:
        regions, ends = cluster_umis_on_position(args.bam_file, args.position_threshold, 
                                             group_method, args.bed_file)
    
    #print(regions)
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
    if group_method=='fromTag':
        bamfilelist = cluster_umis_all_regions(regions, ends, edit_distance_threshold,
                                           args.sample_name, args.bam_file, args.output_path,
                                           args.include_singletons, fasta, bedregions,
                                           num_cpus, args.indel_frequency_threshold,
                                           args.consensus_frequency_threshold,
                                           args.regions_from_tag, starts)
    else:
        bamfilelist = cluster_umis_all_regions(regions, ends, edit_distance_threshold, 
                                           args.sample_name, args.bam_file, args.output_path, 
                                           args.include_singletons, fasta, bedregions, 
                                           num_cpus, args.indel_frequency_threshold, 
                                           args.consensus_frequency_threshold)
    merge_bams(args.output_path, bamfilelist, args.sample_name)
    index_bam_file(args.output_path + '/' + args.sample_name + '_consensus_reads.bam',
              num_cpus)
    consfilelist = [x.rstrip('.bam') + '.cons' for x in bamfilelist]
    merge_cons(args.output_path, consfilelist, args.sample_name)
    cons_file = args.output_path + '/' + args.sample_name + '_cons.tsv'
    if args.remove_large_files:
        os.remove(args.output_path+'/' +args.bam_file)

    statfilelist = [x.rstrip('.bam') + '.hist' for x in bamfilelist]
    merge_stat(args.output_path, statfilelist, args.sample_name)
    duppos = check_duplicate_positions(cons_file)
    if any(duppos):
        merge_duplicate_positions_all_chromosomes(duppos,cons_file,num_cpus)
    merge_duplicate_stat(args.output_path,args.sample_name)
    logging.info("Consensus generation complete, output written to {}, {}".format(args.output_path + 
                 '/' + args.sample_name + '_consensus_reads.bam',
                 args.output_path + '/' + args.sample_name + '_cons.tsv'))

def main(args):
    run_umi_errorcorrect(args)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
