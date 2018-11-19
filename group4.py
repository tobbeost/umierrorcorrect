#!/usr/bin/env python3
import sys
import pysam
from collections import Counter

class Region:
    def __init__(self,pos):
        self.start=pos
        self.end=pos
    def is_inside(self,pos,pos_threshold):
        if pos > self.start - pos_threshold and pos < self.end + pos_threshold:
            if pos < self.start:
                self.start = pos
            if pos > self.end:
                self.end = pos
            return(True)
        else:
            return(False)

def get_chromosome_list_from_bam(f):
    contiglist=[]
    for chrx in f.get_index_statistics():
        if chrx.total > 0:
            contiglist.append(chrx.contig)
    return contiglist

        

def group_by_position(f,chrx,pos_threshold):
    regions={}
    reads=f.fetch(chrx)
    current_pos=-pos_threshold
    current_end=-pos_threshold
    for line in reads:
        pos = line.pos
        if pos > current_end + pos_threshold:
            #new region
            current_pos = pos
            current_end = pos
            regions[pos]=Counter()
            barcode = line.qname.split(':')[-1]
            #if barcode not in regions[current_pos]:
            #    regions[current_pos][barcode]=0
            regions[current_pos][barcode]+=1
        else:
            barcode=line.qname.split(':')[-1]
            #if barcode not in regions[current_pos]:
            #    regions[current_pos][barcode]=0
            regions[current_pos][barcode]+=1
            current_end=pos
    return(regions)
        #if len(regions)==0:
        #    r=Region(pos)
        #    regions.append(r)
        #else:
        #    for rr in regions:
        #        if not rr.is_inside(pos,10):
        #            #new region
        #            rnew=Region(pos)
               
    #for rr in regions:
    #    print(rr.start,rr.end)
def get_max_number_of_barcodes(regions,pos):
    umi,count=regions[pos].most_common(1)[0]
    return(count)

def remove_singleton_regions(regions):
    newregions={}
    for chrx in regions:
        newregions[chrx]=dict((x,y) for (x,y) in regions[chrx].items() if get_max_number_of_barcodes(regions[chrx], x)>1)
    return(newregions)

def readBam(infile):
    '''Read grouped BAM-file (UMI-groups-sorted) and extract sequences from each UMI-group and save one representative read from each group.'''
    with pysam.AlignmentFile(infile,'rb') as f:
        chrs=get_chromosome_list_from_bam(f)
        chrregions={}
        for chrx in chrs:
            chrregions[chrx]=group_by_position(f,chrx,10)

        regions=remove_singleton_regions(chrregions)
        #for chrx in chrs:
        #    print(chrx,regions[chrx])

  
        #for chrx in chrs:
        #    regions2=regions[chrx]
        #    for rr in regions2:
        #        print(chrx,rr,regions2[rr].most_common(10))
    return(regions)

def main(filename):
    regions=readBam(filename)
    for chrx in regions:
        regions2=regions[chrx]
        for rr in regions2:
            print(chrx,rr,regions2[rr].most_common(10))
if __name__=='__main__':
    main(sys.argv[1])
