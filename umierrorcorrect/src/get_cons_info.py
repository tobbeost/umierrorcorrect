#!/usr/bin/env python3

import sys
import pickle
import pysam
# from umi_cluster import umi_cluster
from collections import Counter
from umierrorcorrect.src.get_consensus import get_cons_dict, get_all_consensus, get_reference_sequence
from umierrorcorrect.src.get_regions_from_bed import read_bed, sort_regions, merge_regions, get_annotation, get_annotation2


def get_cons_info(consensus_seq, singletons, fsizes=[0, 1, 2, 3, 4, 5, 7, 10, 20, 30]):
    '''loop through the consensus reads to collapse alleles for each position'''
    cons = {}
    for consensus_read in consensus_seq.values():
        if consensus_read:
            pos = consensus_read.start_pos
            count = consensus_read.count
            if consensus_read.indel_read == 0:
                for base in consensus_read.seq:
                    if pos not in cons:
                        cons[pos] = {}
                    for fsize in fsizes:
                        if fsize == 0:
                            if fsize not in cons[pos]:
                                cons[pos][fsize] = Counter()
                            cons[pos][fsize][base] += count
                        elif count >= fsize:
                            if fsize not in cons[pos]:
                                cons[pos][fsize] = Counter()
                            cons[pos][fsize][base] += 1
                    pos += 1
            else:
                i = 0
                cigar = consensus_read.cigarstring
                for base in consensus_read.seq:
                    c = cigar[i]
                    if pos not in cons:
                        cons[pos] = {}
                    if c == '0':  # match or mismatch
                        for fsize in fsizes:
                            if fsize == 0:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize][base] += count
                            elif count >= fsize:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize][base] += 1
                        pos += 1
                        i += 1
                    elif c == '1':  # insertion
                        for fsize in fsizes:
                            if fsize == 0:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize]['I'] += count
                            elif count >= fsize:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize]['I'] += 1
                        i += 1
                    elif c == '2':  # deletion
                        for fsize in fsizes:
                            if fsize == 0:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize]['D'] += count

                            elif count >= fsize:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize]['D'] += 1
                        deletion = False
                        if cigar[i+1] == '2':
                            deletion = True
                            while deletion:
                                i += 1
                                c = cigar[i]
                                pos += 1
                                if pos not in cons:
                                    cons[pos] = {}
                                for fsize in fsizes:
                                    if fsize == 0:
                                        if fsize not in cons[pos]:
                                            cons[pos][fsize] = Counter()
                                        cons[pos][fsize]['D'] += count

                                    elif count >= fsize:
                                        if fsize not in cons[pos]:
                                            cons[pos][fsize] = Counter()
                                        cons[pos][fsize]['D'] += 1
                                if cigar[i+1] == '2':
                                    deletion = True
                                else:
                                    deletion = False
                        pos += 1
                        if pos not in cons:
                            cons[pos] = {}
                        for fsize in fsizes:
                            if fsize == 0:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize][base] += count

                            elif count >= fsize:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize] = Counter()
                                cons[pos][fsize][base] += 1
                        pos += 1
                        i += 1
    for read in singletons.values():

        if 'I' not in read.cigarstring and 'D' not in read.cigarstring:
            sequence = read.query_sequence
            for qpos, refpos in read.get_aligned_pairs(matches_only=True):
                base = sequence[qpos]
                if refpos not in cons:
                    cons[refpos] = {}
                for fsize in [0, 1]:
                    if fsize not in cons[refpos]:
                        cons[refpos][fsize] = Counter()
                    cons[refpos][fsize][base] += 1

        else:  # insertion or deletion
            positions = read.get_aligned_pairs(matches_only=True)
            q, ref = positions[0]
            q = q - 1
            ref = ref - 1
            for qpos, refpos in positions:
                if not qpos == q + 1:
                    # allele = read.query_sequence[q+1:qpos]
                    # inspos=refpos-1
                    if refpos not in cons:
                        cons[refpos] = {}
                    for fsize in [0, 1]:
                        if fsize not in cons[refpos]:
                            cons[refpos][fsize] = Counter()
                        cons[refpos][fsize]['I'] += 1
                        cons[refpos][fsize][base] += 1
                elif not refpos == ref + 1:
                        # deletion
                    dellength = refpos - (ref + 1)
                    delpos = refpos - dellength
                    if delpos not in cons:
                        cons[delpos] = {}
                    for fsize in [0, 1]:
                        if fsize not in cons[delpos]:
                            cons[delpos][fsize] = Counter()
                        cons[delpos][fsize]['D'] += 1
                sequence = read.query_sequence
                base = sequence[qpos]
                if refpos not in cons:
                    cons[refpos] = {}
                for fsize in [0, 1]:
                    if fsize not in cons[refpos]:
                        cons[refpos][fsize] = Counter()
                    cons[refpos][fsize][base] += 1
                else:
                    sequence = read.query_sequence
                    # qual = read.qual
                    base = sequence[qpos]
                    if refpos not in cons:
                        cons[refpos] = {}
                    for fsize in [0, 1]:
                        if fsize not in cons[refpos]:
                            cons[refpos][fsize] = Counter()
                        cons[refpos][fsize][base] += 1
                q = qpos
                ref = refpos

    return(cons)


def calc_major_nonref_allele_frequency(cons, ref, tot):
    comp = cons
    allele = max(comp, key=comp.get)
    count=cons[allele]
    frac = 1.0*(cons[allele]/tot)
    if frac > 0:
        return((allele, frac, count))
    else:
        return(("", 0, 0))


def write_consensus(f, cons, ref_seq, start, contig, annotation, samplename, only_target_regions):
    bases = ['A', 'C', 'G', 'T', 'I', 'D', 'N']
    # print(list(cons.keys())[0],list(cons.keys())[-1],start,len(ref_seq))

    for pos in sorted(cons):
        
        annotation_pos = get_annotation2(annotation, pos + 1)
        if not (annotation_pos == "" and only_target_regions):
            # if len(ref_seq)<(pos-start+1):
            #     print("error",contig,start,ref_seq)
            refbase = ref_seq[pos - start]

            for fsize in cons[pos]:
                line = []
                line.append(samplename)
                line.append(contig)
                line.append(str(pos + 1))
                line.append(annotation_pos)
                line.append(refbase)
                consline = cons[pos][fsize]
                tot = sum(consline.values())
                nonrefcons = {key: consline[key] for key in consline if key != refbase}
                if len(nonrefcons) > 0:
                    mna, freq, count = calc_major_nonref_allele_frequency(nonrefcons, refbase, tot)
                else:
                    mna = ''
                    freq = 0
                    count = 0
                for base in bases:
                    if base in cons[pos][fsize]:
                        line.append(str(cons[pos][fsize][base]))
                    else:
                        line.append(str(0))
                line.append(str(tot))
                line.append(str(fsize))
                line.append(str(count))
                line.append(str(freq))
                line.append(mna)
                f.write('\t'.join(line) + '\n')


def main(bamfilename, bedfilename):
    with open('/home/xsteto/umierrorcorrect/umi.pickle', 'rb') as f:
        umis = pickle.load(f)
    fasta = pysam.FastaFile('/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta')
    ref_seq = get_reference_sequence(fasta, '17', 7577495, 7577800)
    contig = '17'
    start = 7577495
    end = 7577800
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, contig, start, end, True)
    consensus_seq = get_all_consensus(position_matrix, umis, contig)
    cons = get_cons_info(consensus_seq, singleton_matrix)
    regions = read_bed(bedfilename)
    regions = sort_regions(regions)
    regions = merge_regions(regions, 0)
    annotation = regions[contig]
    with open('out/cons.out', 'w') as f:
        write_consensus(f, cons, ref_seq, start, contig, annotation, False)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
