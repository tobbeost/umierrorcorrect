#!/usr/bin/env python3
'''
UMI error correct, run_umierrorcorrect.py - Run the pipeline

:Author: Tobias Osterlund

Purpose
-------

Run the pipeline

'''

import sys
from umierrorcorrect.src.handle_sequences import read_fastq, read_fastq_paired_end
from umierrorcorrect.preprocess import run_preprocessing
from umierrorcorrect.run_mapping import run_mapping
from umierrorcorrect.umi_error_correct import run_umi_errorcorrect
from umierrorcorrect.src.check_args import check_args_fastq, check_args_bam, get_sample_name
from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics
from umierrorcorrect.call_variants9 import run_call_variants
from umierrorcorrect.version import __version__
import argparse
import os
import logging

def parseArgs():
    parser = argparse.ArgumentParser(description="UmiErrorCorrect v. {}. Pipeline for analyzing  barcoded amplicon sequencing data with \
                                                  Unique molecular identifiers (UMI)".format(__version__))
    group1 = parser.add_argument_group('Input and output files options')
    group1.add_argument('-o', '--output_path', dest='output_path', 
                        help='Path to the output directory, required', required=True)
    group1.add_argument('-r1', '--read1', dest='read1', help='Path to first FASTQ file, R1, required', required=True)
    group1.add_argument('-r2', '--read2', dest='read2', help='Path to second FASTQ file, R2 if applicable')
    group1.add_argument('-r', '--reference', dest='reference_file', 
                        help='Path to the reference sequence in Fasta format (indexed), required', required=True)
    group1.add_argument('-bed', '--bed_file', dest='bed_file', 
                        help='Path to a BED file defining the targeted regions, i.e. chromosomal positions. \
                             The Bed file is used for annotation.')
    group1.add_argument('-s', '--sample_name', dest='sample_name', 
                        help='Sample name that will be used as base name for the output files. \
                              If excluded the sample name will be extracted from the fastq files.')

    group2 = parser.add_argument_group('UMI definition options')
    group2.add_argument('-ul', '--umi_length', dest='umi_length', 
                        help='Length of UMI sequence (number of nucleotides). The UMI is assumed to \
                             be located at the start of the read. Required', required=True)
    group2.add_argument('-sl', '--spacer_length', dest='spacer_length', 
                        help='Length of spacer (The number of nucleotides between the UMI and the beginning \
                             of the read). The UMI + spacer will be trimmed off, and the spacer will be \
                             discarded. Default=%(default)s', default='0')
    group2.add_argument('-mode', '--mode', dest='mode',
                        help="Name of library prep, Either 'single' or 'paired', for single end or paired \
                             end data respectively, [default = %(default)s]", default="paired")
    group2.add_argument('-dual', '--dual_index', dest='dual_index', 
                        help='Include this flag if dual indices are used (UMIs both on R1 and R2)', 
                        action='store_true')
    group2.add_argument('-reverse', '--reverse_index', dest='reverse_index', 
                        help="Include this flag if a single index (UMI) is used, but the UMI is located on R2 \
                             (reverse read). Default is UMI on R1.", action='store_true')
    
    group3 = parser.add_argument_group('UMI clustering options')
    group3.add_argument('-regions_from_bed', dest='regions_from_bed', 
                        help='Include this flag if regions used for UMI clustering and variant calling should be \
                              defined from the BED file. Default is to detect the regions automatically from the BAM file. ',
                              action='store_true')
    group3.add_argument('-d', '--edit_distance', dest='edit_distance_threshold', 
                        help="Edit distance threshold for UMI clustering, [default = %(default)s]",
                        default=1)
    group3.add_argument('-p', '--position_threshold', dest='position_threshold', 
                        help='Position threshold for grouping by position [default = %(default)s]',
                        default=10)
    
    group4 = parser.add_argument_group('Consensus options')
    group4.add_argument('-cons_freq', '--consensus_frequency_threshold', dest='consensus_frequency_threshold',
                        help='Minimum percent of the majority base at a position for consensus to be called. \
                              [default = %(default)s]', default=60.0)
    group4.add_argument('-indel_freq', '--indel_frequency_threshold', dest='indel_frequency_threshold',
                        help='Percent threshold for indels to be included in the consensus read. \
                              [default = %(default)s]', default=60.0)
    group4.add_argument('-singletons', '--include_singletons', dest='include_singletons', action='store_true',
                        help='Include this flag if singleton reads should be included in the output consensus \
                              read bam file. Note that the singletons will not be error corrected')
    
    group5 = parser.add_argument_group('Variant calling options')
    group5.add_argument('-fs','--fsize', dest='fsize', 
                        help='Family size cutoff (consensus cutoff) for variant calling. [default = %(default)s]', 
                        default=3)
    group5.add_argument('-method','--vc-method',dest='vc_method',
                        help="Variant calling method, Either 'count' or 'bbmodel'. [default = %(default)s]", 
                        default='count')
    group5.add_argument('-count','--count_cutoff', dest='count_cutoff',
                        help="Consensus read count cutoff (minimum variant allele depth) for calling a variant if method=count [default = %(default)s]", 
                        default=5)
                        
    group5.add_argument('-Q', '--qscore_cutoff', dest='qvalue_threshold',
                        help='Qscore threshold (Minimum variant significance score) for Variant calling, only if method=bbmodel [default = %(default)s]', 
                        default=20)

    group6 = parser.add_argument_group('Running parameters')
    group6.add_argument('-tmpdir', '--tmp_dir', dest='tmpdir',
                        help="temp directory where the temporary files are written and then removed. \
                              Should be the scratch directory on the node. Default is a temp directory \
                              in the output folder.")
    group6.add_argument('-f', '--force', dest='force',action='store_true',
                        help='Include this flag to force output files to be overwritten')
    group6.add_argument('-t', '--num_threads', dest='num_threads', 
                        help='Number of threads to run the program on. Default=%(default)s', default='1')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMI Error Correct')
    return(args)


def main(args):
    if not args.sample_name:
        args.sample_name = get_sample_name(args.read1, args.mode)
    args=check_args_fastq(args)
    #args=check_args_bam(args)
    fastq_files, nseqs = run_preprocessing(args)  # run preprocessing
    print(fastq_files, nseqs)
    bam_file = run_mapping(args.num_threads, args.reference_file, fastq_files, 
                           args.output_path, args.sample_name)  # run mapping
    args.bam_file = bam_file
    print(args.bam_file)
    run_umi_errorcorrect(args)  #run umi errorcorrect
    cons_bam = args.output_path + '/' + args.sample_name + '_consensus_reads.bam'
    stat_filename = args.output_path + '/' + args.sample_name + '.hist'
    run_get_consensus_statistics(args.output_path, cons_bam, stat_filename, args.sample_name)
    args.cons_file=None
    args.params_file=None
    run_call_variants(args)
    logging.info("Finished UMI Error Correct")

if __name__ == '__main__':
    args=parseArgs()
    main(args)
