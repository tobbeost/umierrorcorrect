#!/usr/bin/env python3
import pysam
import sys
from get_cons_dict2 import consensus_read


def get_stat(consensus_filename, stat_filename):
    with open(stat_filename) as f:
        regions=[]
        for line in f:
            line=line.rstrip()
            regionid, pos, singletons, *rest =line.split('\t')
            singletons=int(singletons.split(": ")[-1])
            regions.append((int(regionid),pos,singletons))
    print(regions)
    hist={}
    with pysam.AlignmentFile(consensus_filename,'rb') as f:
        for regionid,pos,singletons in regions: 
            contig=pos.split(':')[0]
            start=int(pos.split(':')[-1].split('-')[0])
            end=int(pos.split(':')[-1].split('-')[1])
            reads=f.fetch(contig,start,end)
            hist[regionid]=[]
            for read in reads:
                count=int(read.qname.split('=')[-1])
                if count > 1:
                    hist[regionid].append(count)
    print(hist)

#def write_stat(f,consensus_seq,singleton_matrix,contig,start,end):
#    for read in consensus_seq:
#        g.write('{}:{}-{}\t{}'.format(contig,start,end,read.count)
#
#        
#    regions,ends=readBam(bamfilename,position_threshold)

if __name__=='__main__':
    get_stat(sys.argv[1],sys.argv[2])


