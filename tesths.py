#!/usr/bin/env python3
from src.handle_sequences import read_fastq_paired_end

def test():
    r1file='/tmp/tobiastmp/E05_GT1/r190128_102251/E05_GT1_DGT_2.fastq'
    r2file='/tmp/tobiastmp/E05_GT1/r190128_102251/E05_GT1_DGT_1.fastq'
    with open(r1file) as f1, open(r2file) as f2:
        for n1,s1,q1,n2,s2,q2 in read_fastq_paired_end(f1,f2):
            print(n1,n2)


if __name__=='__main__':
    test()
