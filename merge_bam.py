#!/usr/bin/env python3
import sys
import pysam

def merge_bams(output_path,bamfilelist):
    print(output_path+'/consensus_reads.bam')

    with pysam.AlignmentFile(bamfilelist[0],'rb') as f,  pysam.AlignmentFile(output_path+'/consensus_reads.bam','wb',template=f) as g:
        
        for line in f:
            g.write(line)
        for filename in bamfilelist[1:]:
            print(filename)
            with pysam.AlignmentFile(filename,'rb') as f1:
                for line in f1:
                    g.write(line)

def main(bamfilelist):
    output_path='outtest'
    merge_bams(output_path,bamfilelist)

if __name__=='__main__':
    main(sys.argv[1:])   

