#!/usr/bin/env python3
import pysam
import sys
import argparse


def parseArgs():
    parser=argparse.ArgumentParser(description="Filter bam file by removing positions without annotation and with a raw sequencing depth lower than threshold.")
    parser.add_argument('-i', '--infile',dest='infile', help='Path to the input file, required', required=True)
    parser.add_argument('-o', '--outfile',dest='outfile', help='Path to the output file, required', required=True)
    parser.add_argument('-c','--consensus_cutoff', dest='consensus_cutoff',help='Consensus depth cutoff, [default = %(default)s]',default='3')
    args=parser.parse_args(sys.argv[1:])
    return(args)


def filter_bam(infilename,outfilename,consensus_cutoff):
    consensus_cutoff=int(consensus_cutoff)
    with pysam.AlignmentFile(infilename,'rb') as f, pysam.AlignmentFile(outfilename,'wb',template=f) as g:
        reads=f.fetch()
        for read in reads:
            size=int(read.qname.split('=')[-1])
            if size >= consensus_cutoff:
                g.write(read)


if __name__ == '__main__':
    args=parseArgs()
    filter_bam(args.infile, args.outfile, args.consensus_cutoff)
