#!/usr/bin/env python3
import pysam
import sys
import argparse
import matplotlib.pyplot as plt
from umierrorcorrect.src.get_regions_from_bed import read_bed, sort_regions, merge_regions, get_annotation

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)")
    parser.add_argument('-o',  '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-c', '--cons_bam', dest='cons_bam_file', help='Path to the consensus BAM-file')
    parser.add_argument('-hist', '--hist_file', dest='hist_file',
                        help='Path to the .hist file')
    args = parser.parse_args(sys.argv[1:])

    return(args)


class region_cons_stat():
    def __init__(self, regionid, pos, name, singletons):
        self.regionid = regionid
        self.pos = pos
        self.name = name
        self.singletons = singletons

    def add_histogram(self, hist, fsizes):
        total_reads = {}
        umis = {}
        total_reads[0] = sum(hist) + self.singletons
        umis[0] = sum(hist) + self.singletons
        for fsize in fsizes:
            tmp=[x for x in hist if x >=fsize]
            if fsize == 1:
                total_reads[fsize] = sum(tmp) + self.singletons
                umis[fsize] = len(tmp) + self.singletons
            else:
                total_reads[fsize] = sum(tmp)
                umis[fsize] = len(tmp)
        self.hist = hist
        self.total_reads = total_reads
        self.umis = umis
        self.fsizes = fsizes

    def write_stats(self):
        lines=[]
        r0 = self.total_reads[0]
        u0 = self.umis[0]
        line = '\t'.join([str(self.regionid), self.pos, self.name,
                        '0', '1.0', str(r0), str(u0)])
        lines.append(line)
        for fsize in self.fsizes:
            line = '\t'.join([str(self.regionid), self.pos, self.name,
                              str(fsize), str(1.0*(self.total_reads[fsize]/r0)), 
                              str(self.total_reads[fsize]), str(self.umis[fsize])])
            lines.append(line)
        return('\n'.join(lines))

def get_stat(consensus_filename, stat_filename):
    with open(stat_filename) as f:
        regions=[]
        for line in f:
            line=line.rstrip()
            regionid, pos, name, cons, singles, *rest =line.split('\t')
            singles=int(singles.split(": ")[-1])
            regionid=int(regionid)
            regions.append((regionid,pos,singles,name))
            
    print(regions)
    hist={}
    with pysam.AlignmentFile(consensus_filename,'rb') as f:
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
    regionstats=[]
    for regionid,pos,singletons,name in regions:
        stat=region_cons_stat(regionid, pos, name, singletons)
        stat.add_histogram(hist[regionid], fsizes)
        #print(stat.write_stats())
        regionstats.append(stat)
    return(regionstats)


def get_overall_statistics(hist,fsizes):
    histall = region_cons_stat('All','all_regions','',0)
    histall.total_reads = {}
    histall.umis = {}
    fsizesnew=fsizes.copy()
    histall.fsizes = fsizes
    fsizesnew.insert(0,0)
    print(fsizesnew)
    for fsize in fsizesnew:
        histall.total_reads[fsize]=0
        histall.umis[fsize]=0
    
    for region in hist:
        for fsize in fsizesnew:    
            histall.total_reads[fsize] += region.total_reads[fsize]
            histall.umis[fsize] += region.umis[fsize]
    #print(histall.write_stats())
    return(histall)


def plot_histogram(hist,plot_filename):
    umisizesall=[]
    for region in hist:
        umisizesall.extend(hist[region])
    num_bins=1000
    umisizesall.sort(reverse=True)
    print(umisizesall[0:100])
    n,bins,patches=plt.hist(umisizesall,num_bins,facecolor='dodgerblue',alpha=0.5)
    plt.xlabel('Barcode family depth')
    plt.ylabel('Frequency')
    plt.title('Histogram of barcode family depth')
    plt.box(False)
    plt.xlim(0,500)
    plt.savefig(plot_filename)


def run_get_consensus_statistics(output_path, consensus_filename, stat_filename):
    hist=get_stat(consensus_filename,stat_filename)
    fsizes=[1,2,3,4,5,7,10,20,30]
    histall = get_overall_statistics(hist,fsizes)
    print(histall.write_stats())
    for stat in hist:
        print(stat.write_stats())


def main(output_path, consensus_filename, stat_filename):
    run_get_consensus_statistics(output_path,  consensus_filename, stat_filename)


if __name__=='__main__':
    args=parseArgs()
    main(args.output_path, args.cons_bam_file, args.hist_file)


