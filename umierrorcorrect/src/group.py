#!/usr/bin/env python3
import sys
import pysam
from collections import Counter
from umierrorcorrect.src.get_regions_from_bed import read_bed, \
     sort_regions, merge_regions, expand_regions_from_bed


# class Region:
#     def __init__(self, pos):
#         self.start = pos
#         self.end = pos
#     def is_inside(self, pos, pos_threshold):
#     if pos > self.start - pos_threshold and pos < self.end + pos_threshold:
#             if pos < self.start:
#                 self.start = pos
#             if pos > self.end:
#                 self.end = pos
#             return(True)
#         else:
#             return(False)

def get_chromosome_list_from_bam(f):
    contiglist = []
    for chrx in f.get_index_statistics():
        if chrx.total > 0:
            contiglist.append(chrx.contig)
    return contiglist


def group_by_position(f, chrx, pos_threshold):
    regions = {}
    ends = {}
    reads = f.fetch(chrx)
    current_pos = -pos_threshold
    current_end = -(2*pos_threshold) + 1
    # current_aend=-pos_threshold
    for line in reads:
        pos = line.pos
        if pos > current_end + pos_threshold:
            # new region
            if pos != -pos_threshold:
                ends[current_pos] = current_end + 1
            current_pos = pos
            current_end = pos
            # current_aend = line.reference_end
            regions[pos] = Counter()
            barcode = line.qname.split(':')[-1]
            # if barcode not in regions[current_pos]:
            #     regions[current_pos][barcode]=0
            regions[current_pos][barcode] += 1
        else:
            barcode = line.qname.split(':')[-1]
            # if barcode not in regions[current_pos]:
            #     regions[current_pos][barcode]=0
            regions[current_pos][barcode] += 1
            current_end = pos
            # aend=line.reference_end
            # if aend > current_aend:
            #     current_aend = aend
    ends[current_pos] = current_end + 1
    return(regions, ends)
    # if len(regions)==0:
    #     r=Region(pos)
    #     regions.append(r)
    # else:
    #     for rr in regions:
    #         if not rr.is_inside(pos,10):
    #             #new region
    #             rnew=Region(pos)
    # for rr in regions:
#     print(rr.start,rr.end)


def count_umis_in_region(f, chrx, pos_start, pos_end):
    region = Counter()
    reads = f.fetch(chrx, pos_start, pos_end)
    for line in reads:
        # pos = line.pos
        barcode = line.qname.split(':')[-1]
        region[barcode] += 1
    return(region)


def get_max_number_of_barcodes(regions, pos):
    if len(regions[pos]) > 0:
        umi, count = regions[pos].most_common(1)[0]
        return(count)
    else:
        return(0)


def remove_singleton_regions(regions, cutoff):
    newregions = {}
    for chrx in regions:
        newregions[chrx] = dict((x, y) for (x, y) in regions[chrx].items()
                                if get_max_number_of_barcodes(regions[chrx],
                                x) > cutoff)
    return(newregions)


def readBam(infile, position_threshold):
    '''Read grouped BAM-file (UMI-groups-sorted) and extract sequences from
    each UMI-group and save one representative read from each group.'''
    with pysam.AlignmentFile(infile, 'rb') as f:
        chrs = get_chromosome_list_from_bam(f)
        chrregions = {}
        chrends = {}
        for chrx in chrs:
            regions, ends = group_by_position(f, chrx, position_threshold)
            chrregions[chrx] = regions
            chrends[chrx] = ends
            # chrregions[chrx]=group_by_position(f,chrx,position_threshold)

        regions = remove_singleton_regions(chrregions, 2)
        # for chrx in chrs:
        #     print(chrx,regions[chrx])

        # for chrx in chrs:
        #     regions2=regions[chrx]
        #     for rr in regions2:
        #         print(chrx,rr,regions2[rr].most_common(10))
    return(regions, chrends)


def read_bam_from_bed(infile, bedfile, position_threshold):
    chrregions = {}
    chrends = {}
    regions = read_bed(bedfile)
    regions = sort_regions(regions)
    regions = merge_regions(regions, position_threshold)
    regions = expand_regions_from_bed(regions, infile)
    newregions = []
    for contig in regions:
        for a, b, c in regions[contig]:
            newregions.append((contig, a, b, c))

    with pysam.AlignmentFile(infile, 'rb') as f:
        chrs = get_chromosome_list_from_bam(f)
        for contig, start, end, name in newregions:
            if contig in chrs:
                if contig not in chrregions:
                    chrregions[contig] = {}
                chrregions[contig][start] = count_umis_in_region(f, contig,
                                                                 start, end)
                if contig not in chrends:
                    chrends[contig] = {}
                chrends[contig][start] = end
        regions = remove_singleton_regions(chrregions, 2)
    return(regions, chrends)


def read_bam_from_tag(infile):
    regions={}
    starts={}
    ends={}
    with pysam.AlignmentFile(infile, 'rb') as f:
        reads=f.fetch()
        for r in reads:
            contig = r.reference_name
            if contig not in regions:
                regions[contig]={}
                starts[contig]={}
                ends[contig]={}
            try:
                utag = r.get_tag('UG')
            except KeyError:
                print('UG tag not present in BAM file')
                utag = ''
            if utag not in regions[contig]:
                regions[contig][utag] = Counter()
                starts[contig][utag] = r.reference_start
                ends[contig][utag] = r.reference_end
            barcode = r.qname.split(':')[-1]
            regions[contig][utag][barcode] += 1
            if r.reference_end > ends[contig][utag]:
                ends[contig][utag]=r.reference_end
    return(regions,starts,ends)

def main(filename, bedfile):
    position_threshold = 10
    group_method = 'automatic'
    # group_method='fromBed'
    if group_method == 'fromBed':
        regions = read_bam_from_bed(filename, bedfile, position_threshold)
    elif group_method == 'fromTag':
        regions = read_bam_from_tag(filename)
    else:
        regions, chrends = readBam(filename, position_threshold)
    for chrx in regions:
        regions2 = regions[chrx]
        for rr in regions2:
            print(chrx, rr, chrends[chrx][rr], regions2[rr].most_common(10))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
