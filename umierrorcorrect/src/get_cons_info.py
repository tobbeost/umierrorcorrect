#!/usr/bin/env python3

import sys
import pickle
import pysam
# from umi_cluster import umi_cluster
from collections import Counter
from umierrorcorrect.src.get_consensus import get_cons_dict, get_all_consensus, get_reference_sequence
from umierrorcorrect.src.get_regions_from_bed import read_bed, sort_regions, merge_regions, get_annotation


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


def calc_major_nonref_allele_frequency(cons, ref):
    tot = sum(cons.values())
    comp = {key: cons[key] for key in cons if key != ref}
    allele = max(comp, key=comp.get)
    frac = 1.0*(cons[allele]/tot)
    if frac > 0:
        return((allele, frac, tot))
    else:
        return(("", 0, tot))


def write_consensus(f, cons, ref_seq, start, contig, annotation, only_target_regions):
    bases = ['A', 'C', 'G', 'T', 'I', 'D', 'N']
    # print(list(cons.keys())[0],list(cons.keys())[-1],start,len(ref_seq))

    for pos in sorted(cons):
        annotation_pos = get_annotation(annotation, pos)
        if not (annotation_pos == "" and only_target_regions):
            # if len(ref_seq)<(pos-start+1):
            #     print("error",contig,start,ref_seq)
            refbase = ref_seq[pos - start]

            for fsize in cons[pos]:
                line = []
                line.append(contig)
                line.append(str(pos))
                line.append(annotation_pos)
                line.append(refbase)
                if len(cons[pos][fsize]) > 1:
                    mna, freq, tot = calc_major_nonref_allele_frequency(cons[pos][fsize], refbase)
                else:
                    mna = ''
                    freq = 0
                    tot = sum(cons[pos][fsize].values())
                for base in bases:
                    if base in cons[pos][fsize]:
                        line.append(str(cons[pos][fsize][base]))
                    else:
                        line.append(str(0))
                line.append(str(tot))
                line.append(str(fsize))
                line.append(str(freq))
                line.append(mna)
                f.write('\t'.join(line) + '\n')


# def get_all_consensus(position_matrix, umis, contig):
#     consensuses={}
#     for umi in position_matrix:
#         consensuses[umi] = getConsensus3(position_matrix[umi], contig, 50.0, umis[umi])
#     return(consensuses)


# def get_cons_dict(bamfilename, umis, contig, start, end, include_singletons):
#     # print('{}:{}-{}'.format(contig,start,end))
#     position_matrix={}
#     singleton_matrix={}
#     with pysam.AlignmentFile(bamfilename,'rb') as f:
#         alignment=f.fetch(contig,start,end)
#         for read in alignment:
#             barcode=read.qname.split(':')[-1]
#             pos=read.pos
#             if pos>=start and pos <=end:
#             #if barcode not in ['ATGGGGGGGGGG','TCGGGGGGGGGG','TAGGGGGTGGGG','TAGGGGGGGGGG']:
#             #    print('{}:{}-{}'.format(contig,start,end))
#             #    for u in umis:
#             #        print(u, umis[u].centroid, umis[u].count)
#                 cluster=umis[barcode].centroid
#                 cluster_size=umis[barcode].count
#                 if cluster_size>1:
#                     if cluster not in position_matrix:
#                         position_matrix[cluster]=[]
#                     position_matrix[cluster].append(read)
#                 elif include_singletons:
#                     if cluster not in singleton_matrix:
#                         singleton_matrix[cluster]=read
#     return(position_matrix,singleton_matrix)
#
# def write_singleton_reads(singleton_matrix,contig,g):
#     for umi,read in singleton_matrix.items():
#         read.query_name='Singleton_read_{}_{}_Count=1'.format(contig,umi)
#         g.write(read)

def main(bamfilename, bedfilename):
    with open('/home/xsteto/umierrorcorrect/umi.pickle', 'rb') as f:
        umis = pickle.load(f)
    # family_sizes = [0,1,2,3,4,5,7,10,20,30]
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
    # with pysam.AlignmentFile(bamfilename,'rb') as f, pysam.AlignmentFile('consensus_out.bam','wb',template=f) as g:

    #    for cons_read in consensus_seq.values():
    #        if cons_read:
    #            cons_read.write_to_bam(g)
    #    write_singleton_reads(singleton_matrix,'17',g)

    # print(seqs)
    # consensus_seq=get_consensus_seq2(position_matrix)
    # print(position_matrix['CATGGCGAGCAT'])
    # print(umis['CATGGCGAGCCT'].count)

    #    if not consensus_seq[k]:
    #        print(None)
    #    else:
    #        print('@'+k+'_'+str(umis[k].count)+'_'+str(consensus_seq[k].start_pos)+'\n'+consensus_seq[k].seq+'\n+\n'+consensus_seq[k].qual)
    # fasta.close()
    # for k in singleton_matrix:
    #     print(k, umis[k].count, singleton_matrix[k].alignment.qname)

    # print(position_matrix['GCACCCGCGCAC'])
    # print(position_matrix['GCACCCGGGCCC'])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
