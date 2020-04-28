#!/usr/bin/env python
from umierrorcorrect.src.group import readBam, read_bam_from_bed
from umierrorcorrect.src.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from umierrorcorrect.src.get_consensus3 import get_cons_dict, get_all_consensus, write_singleton_reads, get_reference_sequence, get_phred, consensus_read,calc_consensus_probabilities,get_ascii, get_most_common_allele, getConsensus3
from umierrorcorrect.src.get_cons_info import get_cons_info, write_consensus, calc_major_nonref_allele_frequency
from umierrorcorrect.src.get_regions_from_bed import read_bed, sort_regions, merge_regions, get_overlap
from umierrorcorrect.umi_error_correct import cluster_umis_on_position, run_umi_errorcorrect, cluster_consensus_worker, update_bam_header
from umierrorcorrect.version import __version__
import sys
import os
import pysam
from multiprocessing import Pool, cpu_count
import subprocess
import argparse
import logging
import pickle
from collections import Counter


def test_part1(bam_file,reference_fasta):
    contig='chr4'
    pos=54727404
    umi_dict=Counter({'CGGCTCTTCGGT': 189, 'CGGATCTTCGGT': 2, 'CGGCTTTTCGGT': 1, 'CGGCTCGTCGGT': 1})
    ends={'chr4': {-20: -38, 54727404: 54727428}}
    annotations=[(54727426, 54727464, 'KIT_2_v4')]
    argvec=[]
    bamfilelist=[]
    test='CGGCTCTTCGGT'
    i=1
    tmpfilename = 'tmp_{}.bam'.format(i)
    argvec.append((umi_dict, 'test', tmpfilename, int(i), contig, int(pos),
                   int(ends[contig][pos]), int(1), bam_file,
                   False, annotations, reference_fasta, 60,
                   60))
    bamfilelist.append('{}/tmp_{}.bam'.format('failedbams3', i))
    p = Pool(int(1))
    p.map(cluster_consensus_worker, argvec)
    update_bam_header(tmpfilename,'test')
    umi_dict, samplename, tmpfilename, regionid, contig, start, end, edit_distance_threshold, \
    bamfilename, include_singletons, annotations, fasta, indel_frequency_cutoff, \
    consensus_frequency_cutoff = argvec[0]    
    adj_matrix = cluster_barcodes(umi_dict, 1)
    clusters = get_connected_components(umi_dict, adj_matrix)
    umis = merge_clusters(umi_dict, clusters)
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, contig,
                                                      start, end, True)
    seqs={test:position_matrix[test]}
    consensus = getConsensus3(seqs[test],  contig,1,60.0,umis[test],60.0)
    ci = {test:consensus}
    cons_info=get_cons_info(ci, {} , fsizes=[0, 1, 2, 3, 4, 5, 7, 10, 20, 30])
    print(cons_info)
    with pysam.FastaFile(fasta) as f:
        ref_seq = get_reference_sequence(f,'chr4', 54727404, 54727428)
    #with open('cons.out', 'w') as f:
    #    write_consensus(f, cons_info, ref_seq, 0, contig, annotations, 'GIST_test', False)
    return(consensus,cons_info)


def main(bam_file):
    reference_fasta='/medstore/Illumina_Tobias/hg38/hg38.fa'
    consensus, cons_info=test_part1(bam_file,reference_fasta)
  
    consread=consensus
    print(consread.start_pos,consread.splits)
    print(consread.seq)
    print(consread.qual)
    print(consread.cigarstring)
    print(cons_info[54727427][3])
    

if __name__=='__main__':
    bamfile=sys.argv[1]
    main(bamfile)
