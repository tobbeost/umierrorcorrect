#!/usr/bin/env python3

import sys
import pickle
import pysam
from testgroup3 import umi_cluster

def get_consensus(bamfilename,umis):
    with pysam.AlignmentFile(bamfilename,'rb') as f:
        alignment=f.pileup('17',7577495,7577595)
        for pileupcolumn in alignment:
            print(pileupcolumn.pos)



def main(bamfilename):
    with open('/home/xsteto/tmp/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    get_consensus(bamfilename,umis)

if __name__=='__main__':
    main(sys.argv[1])
