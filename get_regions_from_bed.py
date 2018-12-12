#!/usr/bin/env python3
import sys
def read_bed(bedfile):
    regions=[]
    with open(bedfile,'r') as f:
        for line in f:
            contig, start, end, name, *rest = line.split()
            regions.append((contig,int(start),int(end),name))
    return(regions)

def sort_regions(regions):
    newregions=sorted(regions, key=lambda tup: (tup[0],tup[1]))
    return(newregions)
#def get_regions_from_bed(bedfile):

def merge_regions(regions,pos_threshold):
    newregions=[]
    contigs,starts,ends,names=zip(*regions)
    current_start = starts[0]
    current_end = ends[0]
    current_contig = contigs[0]
    current_name = []
    current_name.append(names[0])

    for current_index, start in enumerate(starts):
        if current_end + pos_threshold < start or current_contig != contigs[current_index]:
            newregions.append((current_contig, current_start - pos_threshold, current_end + pos_threshold, ','.join(current_name)))
            current_start = start
            current_contig=contigs[current_index]
            current_name=[]
        current_end=ends[current_index]
        if names[current_index] not in current_name:
            current_name.append(names[current_index])

    newregions.append((current_contig, current_start - pos_threshold, current_end + pos_threshold, ','.join(current_name))) #save the last entry
    return(newregions)

def main(bedfile):
    regions=read_bed(bedfile)
    regions=sort_regions(regions)
    regions=merge_regions(regions,100)
    for r in regions:
        r2=tuple(str(x) for x in r)
        print('\t'.join(r2))

if __name__=='__main__':
    main(sys.argv[1])
