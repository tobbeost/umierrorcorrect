#!/usr/bin/env python
from src.group import readBam, read_bam_from_bed
from src.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from src.get_consensus import get_cons_dict, get_all_consensus, write_singleton_reads, get_reference_sequence
from src.get_cons_info import get_cons_info, write_consensus
from src.get_regions_from_bed import read_bed, sort_regions, merge_regions
import sys
import os
import pysam
from multiprocessing import Pool, cpu_count
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description="Pipeline for analyzing  barcoded amplicon sequencing data with Unique molecular identifiers (UMI)")
    parser.add_argument('-o',  '--output_path', dest='output_path', help='Path to the output directory, required', required=True)
    parser.add_argument('-b', '--bam', dest='bam_file', help='Path to BAM-file')
    parser.add_argument('-bed', '--bed-file', dest='bed_file', help='Path to a BED file defining the targeted regions, i.e. chromosomal positions. The Bed file is used for annotation.')
    parser.add_argument('-regions_from_bed', dest='regions_from_bed', help='Include this flag if regions used for UMI clustering and variant calling should be defined from the BED file. Default is to detect the regions automatically from the BAM file. ', action='store_true')
    parser.add_argument('-r', '--reference', dest='reference_file', help='Path to the reference sequence in Fasta format (indexed), Used for annotation')
    parser.add_argument('-s', '--sample_name', dest='sample_name', help='Sample name that will be used as base name for the output files. If excluded the sample name will be extracted from the BAM file.')
    parser.add_argument('-d', '--edit_distance', dest='edit_distance_threshold', help="Edit distance threshold for UMI clustering, [default = %(default)s]", default=1)
    parser.add_argument('-p', '--position_threshold', dest='position_threshold', help='Position threshold for grouping by position [default = %(default)s]', default=10)
    parser.add_argument('-cons_freq', '--consensus_frequency_threshold', dest='consensus_frequency_threshold', help='Minimum percent of the majority base at a position for consensus to be called. [default = %(default)s]', default=75.0)
    parser.add_argument('-indel_freq', '--indel_frequency_threshold', dest='indel_frequency_threshold', help='Percent threshold for indels to be included in the consensus read. [default = %(default)s]', default=75.0)
    parser.add_argument('-singletons', '--include_singletons', dest='include_singletons', action='store_true', help='Include this flag if singleton reads should be included in the output consensus read bam file. Note that the singletons will not be error corrected')
    parser.add_argument('-t', '--num_threads', dest='num_threads', help='Number of threads to run the program on. If excluded, the number of cpus are automatically detected')
    args = parser.parse_args(sys.argv[1:])
    return(args)
# def umi_cluster_worker(args):
#     umi_dict, region_id, edit_distance_threshold=args
#     adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
#     clusters=get_connected_components(umi_dict,adj_matrix)
#     umis=merge_clusters(umi_dict,clusters)
#     return(umis)
#
#     '''Cluster UMIs based on UMI sequence for one region'''


def check_output_directory(outdir):
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
    umi_dict, tmpfilename, regionid, contig, start, end,  edit_distance_threshold, bamfilename, include_singletons, annotations, fasta, indel_frequency_cutoff = args
    adj_matrix = cluster_barcodes(umi_dict, edit_distance_threshold)
    clusters = get_connected_components(umi_dict, adj_matrix)
    umis = merge_clusters(umi_dict, clusters)
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, contig, start, end, True)  # include_singletons=True
    consensus_seq = get_all_consensus(position_matrix, umis, contig, regionid, indel_frequency_cutoff)
    outfilename = tmpfilename
    with pysam.AlignmentFile(bamfilename, 'rb') as f, pysam.AlignmentFile(outfilename, 'wb', template=f) as g:
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
        if include_singletons:
            write_singleton_reads(singleton_matrix, contig, g)
    cons = get_cons_info(consensus_seq, singleton_matrix)
    consfilename = outfilename.rstrip('.bam') + '.cons'
    statfilename = outfilename.rstrip('.bam') + '.hist'
    if len(cons) > 0:
        startpos = min(list(cons.keys()))
        endpos = max(list(cons.keys())) + 1
        with pysam.FastaFile(fasta) as f:
            ref_seq = get_reference_sequence(f, contig, startpos, endpos)
        with open(consfilename, 'w') as g:
            write_consensus(g, cons, ref_seq, startpos, contig, annotations, False)
    else:  # empty file
        g = open(consfilename, 'w')
        g.close()
    with open(statfilename, 'w') as g2:
        regionname = '{}:{}-{}'.format(contig, start, end)
        g2.write('\t'.join([str(regionid), regionname, 'singletons: ' + str(len(singleton_matrix))]) + '\n')


# def consensus_worker():
# def merge_bams(output_path,bamfilelist):
#    print(output_path+'/consensus_reads.bam')
#    with open(output_path+'/consensus_reads.bam','wb') as g:
#        filename=bamfilelist[0]
#        print(filename)
#        with open(filename,'rb') as f:
#            for line in f:
#                g.write(line)
#        for filename in bamfilelist[1:]:
#            print(filename)
#            with open(filename,'rb') as f:
#                for line in f:
#                    if not line.startswith(b'@'):
#                    #print(line)
#                        g.write(line)

def merge_bams(output_path, bamfilelist, sample_name):
    with pysam.AlignmentFile(bamfilelist[0], 'rb') as f,  pysam.AlignmentFile(output_path + '/' + sample_name + '_consensus_reads.bam', 'wb', template=f) as g:
        for line in f:
            g.write(line)
        if len(bamfilelist) > 1:
            for filename in bamfilelist[1:]:
                # print(filename)
                with pysam.AlignmentFile(filename, 'rb') as f1:
                    for line in f1:
                        g.write(line)
    for filename in bamfilelist:
        os.remove(filename)


def merge_cons(output_path, consfilelist, sample_name):
    with open(output_path + '/' + sample_name + '.cons', 'w') as g:
        g.write('Contig\tPosition\tReference\tA\tC\tG\tT\tI\tD\tN\tCoverage\tConsensus group size\tMax Non-ref Allele Frequency\tMax Non-ref Allele\n')
        for filename in consfilelist:
            with open(filename) as f:
                for line in f:
                    g.write(line)

    for filename in consfilelist:
        os.remove(filename)


def merge_stat(output_path, statfilelist, sample_name):
    with open(output_path + '/' + sample_name + '.hist', 'w') as g:
        for filename in statfilelist:
            with open(filename) as f:
                for line in f:
                    g.write(line)

    for filename in statfilelist:
        os.remove(filename)


def cluster_umis_all_regions(regions, ends, edit_distance_threshold, bamfilename, output_path, include_singletons, fasta, bedregions, num_cpus, indel_frequency_cutoff):
    argvec = []
    bamfilelist = []
    i = 0
    for contig in regions:
        for pos in regions[contig]:
            if contig in bedregions:
                annotations = bedregions[contig]
            else:
                annotations = []
            # endpos = list(regions[contig][pos].keys())[-1]
            tmpfilename = '{}/tmp_{}.bam'.format(output_path, i)
            argvec.append((regions[contig][pos], tmpfilename, int(i), contig, int(pos), int(ends[contig][pos]), int(edit_distance_threshold), bamfilename, include_singletons, annotations, fasta, int(indel_frequency_cutoff)))
            bamfilelist.append('{}/tmp_{}.bam'.format(output_path, i))
            i += 1

    p = Pool(int(num_cpus))
    # print(argvec)
    p.map(cluster_consensus_worker, argvec)
    return(bamfilelist)


def cluster_umis_on_position(bamfilename, position_threshold, group_method, bedfilename=None):
    position_threshold = 20
    # group_method='fromBed'
    # group_method='automatic'
    if group_method == 'fromBed':
        regions, ends = read_bam_from_bed(bamfilename, bedfilename, position_threshold)
    else:
        regions, ends = readBam(bamfilename, position_threshold)
    # for chrx in regions:
    #    regions2=regions[chrx]
    #    for rr in regions2:
    #        print(chrx,rr,regions2[rr].most_common(10))
    return(regions, ends)


def run_umi_errorcorrect(args):
    args.output_path = check_output_directory(args.output_path)
    if args.regions_from_bed:
        group_method = 'fromBed'
    else:
        group_method = 'automatic'
    print(group_method)
    if not args.sample_name:
        args.sample_name = get_sample_name(args.bam_file)
    print(args.bam_file)
    print(args.sample_name)
    regions, ends = cluster_umis_on_position(args.bam_file, args.position_threshold, group_method, args.bed_file)
    nregions = 0
    for chrx in regions:
        nregions += len(regions[chrx])
    print("Number of regions, ", nregions)
    edit_distance_threshold = args.edit_distance_threshold
    if args.num_threads:
        num_cpus = int(args.num_threads)
    else:
        num_cpus = int(cpu_count())
    print(num_cpus)
    fasta = args.reference_file
    bedregions = read_bed(args.bed_file)
    bedregions = sort_regions(bedregions)
    bedregions = merge_regions(bedregions, 0)
    bamfilelist = cluster_umis_all_regions(regions, ends, edit_distance_threshold, args.bam_file, args.output_path, args.include_singletons, fasta, bedregions, num_cpus, args.indel_frequency_threshold)
    merge_bams(args.output_path, bamfilelist, args.sample_name)
    consfilelist = [x.rstrip('.bam') + '.cons' for x in bamfilelist]
    merge_cons(args.output_path, consfilelist, args.sample_name)
    statfilelist = [x.rstrip('.bam') + '.hist' for x in bamfilelist]
    merge_stat(args.output_path, statfilelist, args.sample_name)

def main(args):
    run_umi_errorcorrect(args)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
