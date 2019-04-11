#!/usr/bin/env python3

import os
import re
import subprocess
import errno
def check_output_directory(outdir):
    '''Check if outdir exists, otherwise create it'''
    if os.path.isdir(outdir):
        return(outdir)
    else:
        os.mkdir(outdir)
        return(outdir)

def get_sample_name(filename, mode):
    '''Get the sample name as the basename of the input files.'''
    if mode == 'single':
        sample_name = filename.split('/')[-1].rstrip('fastq').rstrip('fastq.gz')
    elif mode == 'paired':
        sample_name = filename.split('/')[-1].rstrip('fastq').rstrip('fastq.gz')
        if sample_name.endswith('_001'):
            sample_name = sample_name[:-4]
        if re.match('.*R[1-2]$', sample_name):
            sample_name = sample_name[:-2]
        sample_name = sample_name.rstrip('_')
        if re.search('.*_L00[0-9]$',sample_name):
            sample_name=sample_name[:-5]
    elif mode == 'bam':
        sample_name=filename.split('/')[-1]
        if '.sorted' in sample_name:
            sample_name = sample_name.replace('.sorted','')
        sample_name = sample_name.replace('.bam','')
    return(sample_name)

def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name,'--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

def check_args_fastq(args):
    '''Function for checking arguments'''
    args.output_path = check_output_directory(args.output_path)
    is_pigz=is_tool('pigz')
    is_gzip=is_tool('gzip')
    is_bwa=is_tool('bwa')
    if not is_bwa:
        raise ValueError('Cannot find program "bwa". Please install it and add it to the path.')
    if is_pigz:
        args.gziptool='pigz'
    elif is_gzip:
        args.gziptool='gzip'
    else:
        raise ValueError('Cannot find program "gzip" or "pigz". Install one of them and add to the path.')

    #determine the mode (single or paired)
    if not args.read2:
        args.mode = 'single'
    else:
        args.mode = 'paired'
    if not args.sample_name:
        args.sample_name = get_sample_name(args.read1, args.mode)
    if args.dual_index and not args.mode=='paired':
        raise ValueError("Dual index can only be used when both an R1 and R2 file are supplied, exiting.")
    if args.reverse_index and not args.mode=='paired':
        raise ValueError("Reverse index can only be used when both an R1 and R2 file are supplied, exiting")
    try:
        args.umi_length = int(args.umi_length)
    except ValueError as e:
        raise(e + " Barcode length needs to be an integer")
        sys.exit(1)
    try:
        args.spacer_length = int(args.spacer_length)
    except ValueError as e:
        raise(e + " Spacer length needs to be an integer")
        sys.exit(1)
    #check if fastq files exist
    if not os.path.isfile(args.read1):
        raise ValueError("The file specified as r1 ({}) does not exist.".format(args.read1))
    if args.mode == 'paired':
        if not os.path.isfile(args.read2):
            raise ValueError("The file specified as r2 ({}) does not exist.".format(args.read2))
    #heck if umis_in_header file exists
    if args.mode == 'paired':
        f1file=args.output_path + '/' + args.sample_name + '_R1_umis_in_header.fastq.gz'
        f2file=args.output_path + '/' + args.sample_name + '_R2_umis_in_header.fastq.gz'
        if os.path.isfile(f1file) or os.path.isfile(f2file):
            if not args.force:
                raise ValueError("The file {} already exists. Overwrite it by including --force in the command line".format(f1file))
            else:
                os.remove(f1file)
                os.remove(f2file)
    elif args.mode == 'single':
        f1file=args.output_path + '/' + args.sample_name + '_umis_in_header.fastq.gz'
        print(f1file)
        if os.path.isfile(f1file):
            if not args.force:
                raise ValueError("The file {} already exists. Overwrite it by including --force in the command line".format(f1file))
            else:
                os.remove(f1file)
    return(args)
    

def check_args_bam(args):
    '''Function for checking arguments'''
    args.output_path = check_output_directory(args.output_path)
    basenamefile=args.read1
    if not args.sample_name:
        args.sample_name = get_sample_name(basenamefile, args.mode)
    if args.regions_from_bed and not args.bed_file:
        raise ValueError("To use option regions_from_bed a bedfile needs to be supplied, using -bed option")
    return(args)

if __name__=='__main__':
    is_pigz=is_tool('pigz')
    is_gzip=is_tool('gzip')
    is_bwa=is_tool('bwa')
    print(is_pigz)
    print(is_gzip)
    print(is_bwa)
