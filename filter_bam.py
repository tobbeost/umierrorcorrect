#!/usr/bin/env python3
import pysam
import sys

def filter_bam(infilename,outfilename,consensus_cutoff):
    consensus_cutoff=int(consensus_cutoff)
    with pysam.AlignmentFile(infilename,'rb') as f, pysam.AlignmentFile(outfilename,'wb',template=f) as g:
        reads=f.fetch()
        for read in reads:
            size=int(read.qname.split('=')[-1])
            if size >= consensus_cutoff:
                g.write(read)
if __name__=='__main__':
    filter_bam(sys.argv[1],sys.argv[2],sys.argv[3])
