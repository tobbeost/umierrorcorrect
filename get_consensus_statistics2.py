#!/usr/bin/env python3
import pysam
import sys
from get_cons_dict2 import consensus_read
import matplotlib.pyplot as plt

def get_stat(consensus_filename, stat_filename):
    with open(stat_filename) as f:
        regions=[]
        for line in f:
            line=line.rstrip()
            regionid, pos, singles, *rest =line.split('\t')
            singles=int(singles.split(": ")[-1])
            regionid=int(regionid)
            regions.append((regionid,pos,singles))
            
    print(regions)
    hist={}
    with pysam.AlignmentFile(consensus_filename,'rb') as f:
        contig=pos.split(':')[0]
        start=int(pos.split(':')[-1].split('-')[0])
        end=int(pos.split(':')[-1].split('-')[1])
        reads=f.fetch()
        for read in reads:
            idx=read.qname
            if idx.startswith('Consensus_read'):
                regionid=int(idx.split('_')[2])
                count=int(idx.split('=')[-1])
                if regionid not in hist:
                    hist[regionid]=[]
                hist[regionid].append(count)
    #print(hist)
    fsizes=[1,2,3,4,5,7,10,20,30]
    for regionid,pos,singletons in regions:
        statdict=hist[regionid]
        r0=sum(statdict)+singletons
        u0=len(statdict)+singletons
        m0=max(statdict)
        print(regionid,'0','1.0',r0,u0,m0)
        for fsize in fsizes:
            tmp=[x for x in statdict if x >=fsize]
            r1=sum(tmp)
            u1=len(tmp)
            if fsize==1:
                r1=r1+singletons
                u1=u1+singletons
            f1=1.0*(r1/r0)
            print(regionid,fsize,f1,r1,u1)
    return(hist)
#def write_stat(f,consensus_seq,singleton_matrix,contig,start,end):
#    for read in consensus_seq:
#        g.write('{}:{}-{}\t{}'.format(contig,start,end,read.count)
#
#        
#    regions,ends=readBam(bamfilename,position_threshold)
def plot_histogram(hist,plot_filename):
    umisizesall=[]
    for region in hist:
        umisizesall.extend(hist[region])
    num_bins=500
    umisizesall.sort(reverse=True)
    print(umisizesall[0:100])
    n,bins,patches=plt.hist(umisizesall,num_bins,facecolor='dodgerblue',alpha=0.5)
    plt.xlabel('Barcode family depth')
    plt.ylabel('Frequency')
    plt.title('Histogram of barcode family depth')
    plt.box(False)
    plt.xlim(0,100)
    plt.savefig(plot_filename)

def main(consensus_filename,stat_filename):
    hist=get_stat(consensus_filename,stat_filename)
    plot_filename=consensus_filename[:-4]+'.png'
    plot_histogram(hist,plot_filename)

if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])


