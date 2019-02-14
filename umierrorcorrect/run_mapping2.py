#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os


def parseArgs():
    parser = argparse.ArgumentParser(description="Pipeline for analyzing  barcoded amplicon sequencing data with Unique molecular identifiers (UMI)")
    parser.add_argument('-o', '--output_path', dest='output_path', help='Path to the output directory, required', required=True)
    parser.add_argument('-r1', '--read1', dest='read1', help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2', '--read2', dest='read2', help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-r', '--reference', dest='reference_file', help='Path to the reference sequence in Fasta format (indexed), Used for annotation')
    parser.add_argument('-t', '--num_threads', dest='num_threads', help='Number of threads to run the program on. Default=%(default)s', default='1')
    args = parser.parse_args(sys.argv[1:])
    return(args)


def check_output_directory(outdir):
    if os.path.isdir(outdir):
        return(outdir)
    else:
        os.mkdir(outdir)
        return(outdir)


def run_mapping(num_threads, reference_file, fastq_files, output_path):
    output_file = output_path+'/output'
    if len(fastq_files) == 1:
        bwacommand = ['bwa', 'mem', '-t', num_threads, reference_file, fastq_files[0]]
    if len(fastq_files) == 2:
        bwacommand = ['bwa', 'mem', '-t', num_threads, reference_file, fastq_files[0], fastq_files[1]]
    command = ['samtools', 'view', '-Sb', '-@', num_threads, '-o', output_file + '.bam', '-']
    p1 = subprocess.Popen(bwacommand, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(command, stdin=p1.stdout, stdout=subprocess.PIPE)  # pipe bwa output to samtools view
    p2.communicate()
    p1.stdout.close()
    command = ['samtools', 'sort', '-@',  num_threads, output_file + '.bam', '-o', output_file + '.sorted.bam']
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    p.communicate()
    p.wait()
    command = ['samtools', 'index', output_file + '.sorted.bam']
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    p.communicate()
    p.wait()
    return(output_file + '.sorted.bam')


if __name__ == '__main__':
    args = parseArgs()
    args.output_path = check_output_directory(args.output_path)
    if args.read2 is None:
        fastq_files = [args.read1]
    else:
        fastq_files = [args.read1, args.read2]

    bamfile=run_mapping(args.num_threads, args.reference_file, fastq_files, args.output_path)
    print(bamfile)
