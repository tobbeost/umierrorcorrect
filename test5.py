#!/usr/bin/env python3
import gzip
from handle_sequences2 import read_fastq_paired_end
def main(file1,file2):
    with gzip.open(file1,'rb') as r1file, gzip.open(file2,'rb') as r2file:
        for n1,s1,q1,n2,s2,q2 in read_fastq_paired_end(r1file,r2file):
            print(n1,n2)
            print(s1,s2)

if __name__=='__main__':
    import sys
    main(sys.argv[1],sys.argv[2])
