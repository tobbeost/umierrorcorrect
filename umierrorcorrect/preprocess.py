#!/usr/bin/env python3
'''
UMI error correct, preprocess.py - remove UMI and append to the header of Fastq sequences.
================================

:Author: Tobias Osterlund

Purpose
-------

Preprocess the fastq files by removing the unique molecular index and add it to the header of the fastq entry.

'''
import sys
import gzip
from umierrorcorrect.src.handle_sequences import read_fastq, read_fastq_paired_end
from umierrorcorrect.src.check_args import check_args_fastq
import argparse
import os
import logging
import subprocess


def parseArgs():
    parser = argparse.ArgumentParser(description="UmiErrorCorrect v.{} Pipeline for analyzing  barcoded amplicon sequencing data with \
                                                  Unique molecular identifiers (UMI)")
    parser.add_argument('-o', '--output_path', dest='output_path', help='Path to the output directory, required', 
                         required=True)
    parser.add_argument('-r1', '--read1', dest='read1', help='Path to first FASTQ file, R1, required', 
                         required=True)
    parser.add_argument('-r2', '--read2', dest='read2', help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-ul', '--umi_length', dest='umi_length', 
                        help='Length of UMI sequence (number of nucleotides).  \
                              The UMI is assumed to be located at the start of \
                              the read. Required', required=True)
    parser.add_argument('-sl', '--spacer_length', dest='spacer_length', 
                        help='Length of spacer (The number of nucleotides between the UMI and the beginning of the read). \
                              The UMI + spacer will be trimmed off, and the spacer will be discarded. Default=%(default)s', 
                              default='0')
    #parser.add_argument('-mode', '--mode', dest='mode', 
    #                    help="Name of library prep, Either 'single' or 'paired', for single end or paired end data \
    #                          respectively, [default = %(default)s]", default="paired")
    parser.add_argument('-dual', '--dual_index', dest='dual_index', 
                        help='Include this flag if dual indices are used (UMIs both on R1 and R2)', 
                        action='store_true')
    parser.add_argument('-reverse', '--reverse_index', dest='reverse_index', 
                        help="Include this flag if a single index (UMI) is used, but the UMI is located on R2 \
                              (reverse read). Default is UMI on R1.", action='store_true')
    parser.add_argument('-s', '--sample_name', dest='sample_name', 
                        help='Sample name which is used as base name \
                              for the output files. If excluded the sample name is automatically extracted from the \
                              name of the fastq file(s).')
    parser.add_argument('-tmpdir', '--tmp_dir', dest='tmpdir', 
                        help="temp directory where the temporary files are written and then removed. Should be the \
                              scratch directory on the node. Default is a temp directory in the output folder.")
    parser.add_argument('-f', '--force', dest='force',action='store_true',
                        help='Include this flag to force output files to be overwritten')
    parser.add_argument('-t', '--num_threads', dest='num_threads', 
                        help='Number of threads to run the program on. Default=%(default)s', default='1')
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


def generate_random_dir(tmpdir):
    '''Generate a directory for storing temporary files, using a timestamp.'''
    import datetime
    newtmpdir = tmpdir + '/r' + datetime.datetime.now().strftime("%y%m%d_%H%M%S") + '/'
    newtmpdir = check_output_directory(newtmpdir)
    return(newtmpdir)


def run_unpigz(filename, tmpdir, num_threads, program):
    '''Unzip the fastq.gz files using parallel gzip (pigz).'''
    outfilename = tmpdir + '/' + filename.split('/')[-1].rstrip('.gz')
    if program=='pigz':
        command = ['unpigz', '-p',  num_threads, '-c', filename]
    elif program=='gzip':
        command = ['gunzip', '-c', filename]
    with open(outfilename, 'w') as g:
        p = subprocess.Popen(command, stdout=g)
        p.communicate()
        p.wait()
    return(outfilename)

def run_gunzip(filename, tmpdir):
    '''Unzip the fastq.gz files using parallel gzip (pigz).'''
    outfilename = tmpdir + '/' + filename.split('/')[-1].rstrip('.gz')
    command = ['gunzip', '-c', filename]
    with open(outfilename, 'w') as g:
        p = subprocess.Popen(command, stdout=g)
        p.communicate()
        p.wait()
    return(outfilename)


def run_pigz(filename, num_threads, program):
    '''Zip the new fastq files with parallel gzip (pigz).'''
    if program=='pigz':
        command = ['pigz', '-p', num_threads, filename]
    elif program=='gzip':
        command = ['gzip', filename]
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    p.communicate()
    p.wait()



def preprocess_se(infilename, outfilename, barcode_length, spacer_length):
    '''Run the preprocessing for single end data (one fastq file).'''
    with open(infilename) as f, open(outfilename, 'w') as g:
        read_start = barcode_length + spacer_length
        nseqs = 0
        for name, seq, qual in read_fastq(f):
            nseqs += 1
            barcode = seq[:barcode_length]
            # g.write(name+':'+barcode+'\n'+rest+'\n'+qualname+'\n'+qual[12+11:]+'\n')
            parts = name.split()
            newname = ':'.join([parts[0], barcode]) + ' ' + parts[-1]
            g.write('\n'.join([newname, seq[read_start:], '+', qual[read_start:]]) + '\n')
    return(nseqs)

def preprocess_pe(r1file, r2file, outfile1, outfile2, barcode_length, spacer_length, dual_index):
    '''Run the preprocessing for paired end data (two fastq files).'''
    read_start = barcode_length + spacer_length
    with open(r1file) as f1, open(r2file) as f2, open(outfile1, 'w') as g1, open(outfile2, 'w') as g2:
        nseqs = 0
        for name1, seq1, qual1, name2, seq2, qual2 in read_fastq_paired_end(f1, f2):
            nseqs += 1
            if dual_index:
                barcode = seq1[:barcode_length] + seq2[:barcode_length]
            else:
                barcode = seq1[:barcode_length]
            parts1 = name1.split()
            parts2 = name2.split()
            newname1 = ':'.join([parts1[0], barcode]) + ' ' + parts1[-1]
            newname2 = ':'.join([parts2[0], barcode]) + ' ' + parts2[-1]
            g1.write('\n'.join([newname1, seq1[read_start:], '+', qual1[read_start:]]) + '\n')
            if dual_index:
                g2.write('\n'.join([newname2, seq2[read_start:], '+', qual2[read_start:]]) + '\n')
            else:
                g2.write('\n'.join([newname2, seq2, '+', qual2]) + '\n')
    return(2*nseqs)

def run_preprocessing(args):
    '''Start preprocessing.'''
    logging.info("Start preprocessing of sample {}".format(args.sample_name))

    if args.tmpdir:
        newtmpdir = generate_random_dir(args.tmpdir)
    else:
        newtmpdir = generate_random_dir(args.output_path)
    # args.chunksize=int(args.chunksize)
    # Unzip the fastq.gz files
    if not args.read1.endswith('gz'):
        r1file = args.read1
        removerfiles = False
        if args.mode == 'paired':
            r2file = args.read2
    else:
        removerfiles = True
        if args.mode == 'paired':
            r1file = run_unpigz(args.read1, newtmpdir, args.num_threads, args.gziptool)
            r2file = run_unpigz(args.read2, newtmpdir, args.num_threads, args.gziptool)
        else:
            r1file = run_unpigz(args.read1, newtmpdir, args.num_threads, args.gziptool)

    logging.info('Writing output files to {}'.format(args.output_path))
    if args.mode == 'single':
        outfilename = args.output_path + '/' + args.sample_name + '_umis_in_header.fastq'
        nseqs = preprocess_se(r1file, outfilename, args.umi_length, args.spacer_length)
        run_pigz(outfilename, args.num_threads, args.gziptool)
        os.remove(r1file)
        os.rmdir(newtmpdir)
        fastqfiles=[outfilename + '.gz']
    else:
        if args.reverse_index:
            # switch forward and reverse read
            r1filetmp = r1file
            r1file = r2file
            r2file = r1filetmp
            outfile1 = args.output_path + '/' + args.sample_name + '_R2_umis_in_header.fastq'
            outfile2 = args.output_path + '/' + args.sample_name + '_R1_umis_in_header.fastq'
        else:
            # r1file=args.read1
            # r2file=args.read2
            outfile1 = args.output_path + '/' + args.sample_name + '_R1_umis_in_header.fastq'
            outfile2 = args.output_path + '/'+ args.sample_name + '_R2_umis_in_header.fastq'
        nseqs = preprocess_pe(r1file, r2file, outfile1, outfile2, args.umi_length, args.spacer_length, args.dual_index)
        run_pigz(outfile1, args.num_threads,args.gziptool)
        run_pigz(outfile2, args.num_threads,args.gziptool)
        if removerfiles == True and os.path.isfile(r1file):
            os.remove(r1file)
        if removerfiles == True and os.path.isfile(r2file):
            os.remove(r2file)
        os.rmdir(newtmpdir)
        fastqfiles=[outfile1 + '.gz', outfile2 + '.gz']
    logging.info("Finished preprocessing")
    return(fastqfiles, nseqs)

def main(args):
    try:
        args = check_args_fastq(args)  # check if combination of arguments are correct
    except ValueError as e:
        print(e)
        sys.exit(1)
    fastqfiles, nseqs = run_preprocessing(args)
    print(nseqs)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
