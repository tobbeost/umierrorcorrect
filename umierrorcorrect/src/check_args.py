#!/usr/bin/env python3

import os

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
        sample_name = filename.split('/')[-1].rstrip('fastq').rstrip('fastq.gz').rstrip('_R012')
    elif mode == 'bam':
        sample_name=filename.split('/')[-1]
        if '.sorted' in sample_name:
            sample_name = sample_name.replace('.sorted','')
        sample_name = sample_name.replace('.bam','')
    return(sample_name)


def check_args_fastq(args):
    '''Function for checking arguments'''
    args.output_path = check_output_directory(args.output_path)
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
