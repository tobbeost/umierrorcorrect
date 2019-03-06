#!/usr/bin/env python3
import sys


def read_bed(bedfile):
    regions = {}
    with open(bedfile, 'r') as f:
        for line in f:
            line = line.strip()
            parts = line.split()
            contig, start, end, name  = parts[0:4]
            if contig not in regions:
                regions[contig] = []
            regions[contig].append((int(start), int(end), name))
    return(regions)


def sort_regions(regions):
    newregions = {}
    for contig in regions:
        newregions[contig] = sorted(regions[contig], key=lambda tup: tup[0])
    return(newregions)


def merge_regions(regions, pos_threshold):
    newregions = {}
    for contig in regions:
        newregions[contig] = []
        starts, ends, names = zip(*regions[contig])
        current_start = starts[0]
        current_end = ends[0]
        current_name = []
        current_name.append(names[0])

        for current_index, start in enumerate(starts):
            if current_end + pos_threshold < start:
                newregions[contig].append((current_start - pos_threshold,
                                           current_end + pos_threshold,
                                           ','.join(current_name)))
                current_start = start
                current_name = []
            current_end = ends[current_index]
            if names[current_index] not in current_name:
                current_name.append(names[current_index])

        newregions[contig].append((current_start - pos_threshold,
                                   current_end + pos_threshold,
                                   ','.join(current_name)))  # save last entry
    return(newregions)


def get_annotation(regions, pos):
    for start, end, name in regions:
        if pos >= start and pos <= end:
            return(name)
            break
    else:
        return("")


def get_overlap(annotation_regions, start, end):
    for annotation_start, annotation_end, annotation_name in annotation_regions:
        if annotation_start <= end and start <= annotation_end:  # test for overlap
            return(annotation_name)
            break
    else:
        return("")


def main(bedfile):
    regions = read_bed(bedfile)
    regions = sort_regions(regions)
    regions = merge_regions(regions, 0)
    for c in regions:
        for r in regions[c]:
            r2 = tuple(str(x) for x in r)
            print(c, '\t'.join(r2))
    contig = '17'
    r = regions[contig]
    for pos in range(7577497, 7577600):
        print(pos, get_annotation(r, pos))


if __name__ == '__main__':
    main(sys.argv[1])
