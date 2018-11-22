#!/usr/bin/env python3
import gzip
from handle_sequences2 import read_fastq
import time
def main(filename):
    seqdict={}
    
    start_time = time.time()
    
    with gzip.open(filename,'rb') as f:
        for name,seq,qual in read_fastq(f):
            seqdict[name]=seq

    print("--- %s seconds ---" % (time.time() - start_time))
    print(seqdict[b'@MN00417:24:000H2F35Y:1:11102:7857:1052 1:N:0:8'])
    a=list(seqdict.keys())[0]
    print(a.split())
    #seqdict={}

    #start_time = time.time()
    #for name,seq,qual in readfq(fp):
    #    seqdict[name]=seq

    #print("--- %s seconds ---" % (time.time() - start_time))

if __name__=='__main__':
    import sys
    main(sys.argv[1])
