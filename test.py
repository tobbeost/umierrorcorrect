#!/usr/bin/env python3
import pickle
import pysam
from testgroup3 import umi_cluster


def main():
    filename='/home/xsteto/Cons_depth_estimation/20ng-10-depth-1_S13/output.sorted.bam'
    clustlist=['GCACCCGCGCCC','GCACCCGCGCAC','GCACCCGGGCCC']
    key=clustlist[0]
    with open('/home/xsteto/tmp/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    umig=umis[key]
    print(umig.centroid)
    print(umig.count)
    n=0
    #with pysam.AlignmentFile(filename,'rb') as f:
    #    reads=f.fetch('17',7577497,7577600)
    #    for read in reads:
    #        barcode = read.qname.split(':')[-1]
    #        if barcode in clustlist:
    #            n+=1
    print(n)
    with pysam.AlignmentFile(filename,'rb') as f:
        alignment=f.pileup('17',7577497,7577600, max_depth=1000000)
        for pileupcolumn in alignment:
            for read in pileupcolumn.pileups:
                barcode = read.alignment.qname.split(':')[-1]
                if barcode in clustlist and pileupcolumn.pos==7577500:
                    n+=1
    print(n)
if __name__=='__main__':
    main()
