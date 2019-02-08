#!/usr/bin/env python3

import sys
import pickle
import pysam
# from umi_cluster import umi_cluster
# from collections import Counter
from math import log10
from itertools import groupby


class consensus_read:

    '''Class for representing a consensus read, useful for writing to BAM'''
    def __init__(self, contig, regionid, position_start, name, count):
        self.contig = contig
        self.start_pos = position_start
        self.seq = ''
        self.qual = ''
        self.indel_read = 0
        self.nmtag = 0
        self.cigarstring = ''
        self.name = 'Consensus_read_{}_{}_Count={}'.format(regionid,
                                                           name, count)
        self.count = count

    def add_base(self, base, qual):
        self.seq = self.seq+base
        self.qual = self.qual+qual
        self.cigarstring += '0'  # todo mismatch

    def add_insertion(self, sequence):
        self.seq = self.seq + sequence
        self.qual += ']'*len(sequence)
        self.cigarstring += '1'*len(sequence)
        self.nmtag += len(sequence)
        self.indel_read = 1

    def add_deletion(self, dellength):
        self.cigarstring += '2'*dellength
        self.nmtag += int(dellength)
        self.indel_read = -1

    def get_cigar(self):
        groups = groupby(self.cigarstring)
        cigar = tuple((int(label),
                       sum(1 for _ in group)) for label, group in groups)
        return(cigar)

    def write_to_bam(self, f):
        a = pysam.AlignedSegment()
        a.query_name = self.name
        a.query_sequence = self.seq
        a.flag = 0
        a.reference_id = f.references.index(self.contig)
        a.reference_start = self.start_pos
        a.mapping_quality = 60
        a.cigar = self.get_cigar()
        a.query_qualities = pysam.qualitystring_to_array(self.qual)
        a.tags = (("NM", self.nmtag),
                  ("RG", "L1"))
        f.write(a)


def get_reference_sequence(fasta, chrx, start, stop):
    ref = fasta.fetch(chrx, start, stop).upper()
    return(ref)


def get_most_common_allele(cons_pos):
    cons_dict = {}
    for x, y in cons_pos.items():
        if x in 'ID':
            for al in y:
                cons_dict[x+'{}'.format(al)] = y[al]
        else:
            cons_dict[x] = len(y)
    cons_allele = max(cons_dict, key=cons_dict.get)  # get highest count allele
    cons_percent = (cons_dict[cons_allele] / sum(cons_dict.values())) * 100
    return(cons_allele, cons_percent)


# def get_consensus_seq2(position_matrix):
#     consensus_sequences = {}
#     for barcode in position_matrix:
#         consread = consensus_read('17', list(position_matrix[barcode].keys())[0])
#         for pos in position_matrix[barcode]:
#             if 'I' in position_matrix[barcode][pos] or 'D' in position_matrix[barcode][pos]:
#                 cons_allele, cons_percent = get_most_common_allele(position_matrix[barcode][pos])
#                 if cons_allele in 'ID':
#                     if cons_percent >= 50.0:
#                         consread.add_base(cons_allele, get_ascii(60))
#                         consread.indel(cons_allele)
#                     else:
#                         consread.add_base('N', get_ascii(0))
#                 else:
#                     cons_base,cons_qual = calc_consensus_probabilities(position_matrix[barcode][pos])
#                     consread.add_base(cons_base, get_ascii(cons_qual))
#             else:
#                 cons_base,cons_qual=calc_consensus_probabilities(position_matrix[barcode][pos])
#                 consread.add_base(cons_base,get_ascii(cons_qual))
#         consensus_sequences[barcode]=consread
#     return(consensus_sequences)

# def check_if_match(read,ref):
#     if read==ref:
#         return('{}M'.format(len(read)))
#     else:
#         return('X')
#
# def get_cigar_string(cons,ref):
#     cons_sequence=''.join(cons.values())
#     cons_sequence=cons
#     if not 'D' in cons_sequence and not 'I' in cons_sequence:
#     if 'D' in cons_sequence:
#         parts=cons_sequence.split('D')
#         counter=0
#         dellength=1
#         cigar=''
#         for part in parts:
#             if part not in '':
#                 cigar=cigar+'{}D'.format(dellength)
#                 dellength=1
#                 cigar=cigar+'{}M'.format(len(part))
#                 counter=counter+len(part)+1
#             else:
#                 counter+=1
#                 dellength+=1
#     cigar=cigar[2:]
#     return(cigar)
#
#
#     else:
#         if cons_sequence==ref:
#             return('{}M'.format(len(cons)))
#         else:


def get_phred(character):
    '''Returns the numeric value of ASCII character associated with phred score (offset 33)'''
    value = ord(character) - 33
    return(value)


def get_ascii(value):
    ascii_letter = chr(value + 33)
    return(ascii_letter)


def calc_consensus(base, cons_pos):
    prod = 1
    for nucl in cons_pos:
        if nucl in base:
            for phred in cons_pos[nucl]:
                prod = prod*(1 - (10**(-phred/10)))
        if nucl not in base and nucl in 'ATCG':
            for phred in cons_pos[nucl]:
                prod = prod*(10**(-phred/10))
    return(prod)


def calc_consensus_probabilities(cons_pos):
    p = {base: calc_consensus(base, cons_pos) for base in 'ATCG'}
    denom = sum(p.values())
    if denom > 0:
        probs = {base: p[base]/denom for base in 'ATCG'}
    else:
        probs = {base: 0 for base in 'ATCG'}
    cons_base = max(probs, key=probs.get)
    if probs[cons_base] == 1:
        cons_phred = 60
    else:
        cons_phred = round(-10*log10(1-probs[cons_base]))
        if cons_phred > 60:
            cons_phred = 60
    return(cons_base, cons_phred)


# def get_consensus(bamfilename, umis, family_sizes, ref_seq):
#     consensus_sequence = {}
#     for fs in family_sizes:
#         consensus_sequence[fs] = {}
#     position_matrix = {}
#     singleton_matrix = {}
#     with pysam.AlignmentFile(bamfilename, 'rb') as f:
#         alignment = f.pileup('17', 7577497, 7577600,max_depth=1000000)
#         n=0
#         clustlist=['GCACCCGCGCCC','GCACCCGCGCAC','GCACCCGGGCCC']
#         for pileupcolumn in alignment:
#             pos=pileupcolumn.pos
#             #ref_base=ref_seq[pos-7577497]
#             #print(pos)
#             for read in pileupcolumn.pileups:
#                 barcode=read.alignment.qname.split(':')[-1]
#                 #if barcode in clustlist and pos==7577513:
#                 #    n+=1
#                 cluster=umis[barcode].centroid
#                 cluster_size=umis[barcode].count
#                 if cluster_size > 1:
#                     if not read.is_del and not read.indel:
#                         q=get_phred(read.alignment.qual[read.query_position])
#                         al=read.alignment.query_sequence[read.query_position]
#                         if cluster in 'CATGGCGAGCCT' and pos==7577503:
#                             print(allele)
#                             print(barcode)
#                             print(cluster)
#                             n+=1
#                     elif read.indel>0: #insertion
#                         #print(read.indel)
#                         al='I'+read.alignment.query_sequence[read.query_position:read.query_position+abs(read.indel)+1]
#                         q=get_phred(read.alignment.qual[read.query_position])
#                         #al=read.alignment.query_sequence[read.query_position:read.query_position+abs(read.indel)+1]
#                     elif read.indel<0: #deletion
#                         #print(read.indel)
#                         #al=read.alignment.query_sequence[read.query_position+1]
#                         #ref_base=ref_seq[(pos-7577497):(pos-7577497)+abs(read.indel)]
#                         al='D'+str(abs(read.indel))
#                         q=get_phred(read.alignment.qual[read.query_position])
#                     if cluster not in position_matrix:
#                         position_matrix[cluster]={}
#                     if pos not in position_matrix[cluster]:
#                         position_matrix[cluster][pos]={}
#                     #allele=(ref_base,al)
#                     allele=al
#                     if allele not in position_matrix[cluster][pos]:
#                         position_matrix[cluster][pos][allele]=[]
#                     position_matrix[cluster][pos][allele].append(q)
#                 else:
#                     if cluster not in singleton_matrix:
#                         singleton_matrix[cluster]=read
#    print(n)
#    return(position_matrix,singleton_matrix)


def getConsensus3(group_seqs, contig, regionid, indel_freq_threshold, umi_info):
    '''Takes a list of pysam entries (rows in the BAM file) as input and generates a consensus sequence.'''
    consensus = {}
    for read in group_seqs:
        if read.cigarstring:
            if 'I' not in read.cigarstring and 'D' not in read.cigarstring:  # no indel in read
                sequence = read.seq
                qual = read.qual
                for qpos, refpos in read.get_aligned_pairs(matches_only=True):
                    base = sequence[qpos]
                    if refpos not in consensus:
                        consensus[refpos] = {}
                    if base not in consensus[refpos]:
                        consensus[refpos][base] = []
                    consensus[refpos][base].append(get_phred(qual[qpos]))
            else:  # indel at next position
                positions = read.get_aligned_pairs(matches_only=True)
                q, ref = positions[0]
                q = q - 1
                ref = ref - 1
                for qpos, refpos in positions:
                    if not qpos == q+1:
                        # insertion
                        allele = read.seq[q+1:qpos]
                        # inspos=refpos-1
                        if refpos not in consensus:
                            consensus[refpos] = {}
                        if 'I' not in consensus[refpos]:
                            consensus[refpos]['I'] = {}
                        if allele not in consensus[refpos]['I']:
                            consensus[refpos]['I'][allele] = 0
                        consensus[refpos]['I'][allele] += 1
                        sequence = read.seq
                        qual = read.qual
                        base = sequence[qpos]
                        if base not in consensus[refpos]:
                            consensus[refpos][base] = []
                        consensus[refpos][base].append(get_phred(qual[qpos]))
                    elif not refpos == ref + 1:
                        # deletion
                        dellength = refpos - (ref+1)
                        delpos = refpos - dellength
                        if delpos not in consensus:
                            consensus[delpos] = {}
                        if 'D' not in consensus[delpos]:
                            consensus[delpos]['D'] = {}
                        if dellength not in consensus[delpos]['D']:
                            consensus[delpos]['D'][dellength] = 0
                        consensus[delpos]['D'][dellength] += 1
                        sequence = read.seq
                        qual = read.qual
                        base = sequence[qpos]
                        if refpos not in consensus:
                            consensus[refpos] = {}
                        if base not in consensus[refpos]:
                            consensus[refpos][base] = []
                        consensus[refpos][base].append(get_phred(qual[qpos]))
                    else:
                        sequence = read.seq
                        qual = read.qual
                        base = sequence[qpos]
                        if refpos not in consensus:
                            consensus[refpos] = {}
                        if base not in consensus[refpos]:
                            consensus[refpos][base] = []
                        consensus[refpos][base].append(get_phred(qual[qpos]))
                    q = qpos
                    ref = refpos
            # curr_index=0
            # for typex,num_bases in read.cigar:
            #     if typex==0:
            #         for i in range(curr_index,curr_index+num_bases):
            #             print(readdict[i])
            # print(read.to_string())#consensusseq={}
    if len(consensus) > 0:
        consensus_sorted = sorted(consensus)
        consread = consensus_read(contig, regionid, consensus_sorted[0], umi_info.centroid, umi_info.count)
        for pos in sorted(consensus_sorted):
            if 'I' in consensus[pos]:
                # first add the insertion if it is in the majority of the reads, then add the base at the next position
                cons_dict = consensus[pos]['I']
                cons_allele = max(cons_dict, key=cons_dict.get)
                cons_num = cons_dict[cons_allele]
                percent = (cons_num / len(group_seqs))*100.0
                if percent >= indel_freq_threshold:
                    sequence = cons_allele
                    consread.add_insertion(sequence)
                del(consensus[pos]['I'])
                cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                consread.add_base(cons_base, get_ascii(cons_qual))

            elif 'D' in consensus[pos]:
                # add the deletions
                a, percent = get_most_common_allele(consensus[pos])
                if a.startswith('D'):
                    if percent >= indel_freq_threshold:
                        dellength = int(a.lstrip('D'))
                        consread.add_deletion(dellength)
                    else:
                        consread.add_base('N', get_ascii(0))
                elif percent >= indel_freq_threshold:
                    cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                    consread.add_base(cons_base, get_ascii(cons_qual))
                else:
                    consread.add_base('N', get_ascii(0))
            else:
                cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                consread.add_base(cons_base, get_ascii(cons_qual))

        return(consread)
    else:
        return(None)

    #    for position in consensus:
    #     if len(consensus[position])==1: #only one option
    #         consensusseq[position]=list(consensus[position].keys())[0]
    #     else:
    #         cons_allele = max(consensus[position], key = consensus[position].get)
    #         cons_denom = sum(consensus[position].values())
    #         cons_percent = (consensus[position][cons_allele]/sum(consensus[position].values())) * 100
    #         if cons_percent>=freq_cutoff:
    #             consensusseq[position]=cons_allele
    #         else:
    #             consensusseq=None
    #             break
    # return(consensusseq)


def get_all_consensus(position_matrix, umis, contig, regionid, indel_frequency_cutoff):
    consensuses = {}
    for umi in position_matrix:
        consensuses[umi] = getConsensus3(position_matrix[umi], contig, regionid, indel_frequency_cutoff, umis[umi])
    return(consensuses)


def get_cons_dict(bamfilename, umis, contig, start, end, include_singletons):
    # print('{}:{}-{}'.format(contig,start,end))
    position_matrix = {}
    singleton_matrix = {}
    with pysam.AlignmentFile(bamfilename, 'rb') as f:
        alignment = f.fetch(contig, start, end)
        for read in alignment:
            barcode = read.qname.split(':')[-1]
            pos = read.pos
            if pos >= start and pos <= end:
                # if barcode not in ['ATGGGGGGGGGG','TCGGGGGGGGGG','TAGGGGGTGGGG','TAGGGGGGGGGG']:
                #     print('{}:{}-{}'.format(contig,start,end))
                #     for u in umis:
                #         print(u, umis[u].centroid, umis[u].count)
                cluster = umis[barcode].centroid
                cluster_size = umis[barcode].count
                if cluster_size > 1:
                    if cluster not in position_matrix:
                        position_matrix[cluster] = []
                    position_matrix[cluster].append(read)
                elif include_singletons:
                    if cluster not in singleton_matrix:
                        singleton_matrix[cluster] = read
    return(position_matrix, singleton_matrix)


def write_singleton_reads(singleton_matrix, contig, g):
    for umi, read in singleton_matrix.items():
        read.query_name = 'Singleton_read_{}_{}_Count=1'.format(contig, umi)
        g.write(read)


def main(bamfilename):
    with open('/home/xsteto/umierrorcorrect/umi.pickle', 'rb') as f:
        umis = pickle.load(f)
    # family_sizes=[0,1,2,3,4,5,7,10,20,30]
    # fasta=pysam.FastaFile('/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta')
    # ref_seq=get_reference_sequence(fasta,'17',7577497,7577800)
    contig = '17'
    start = 7577497
    end = 7577800
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, contig, start, end)
    consensus_seq = get_all_consensus(position_matrix, umis, contig)
    with pysam.AlignmentFile(bamfilename, 'rb') as f, pysam.AlignmentFile('consensus_out.bam', 'wb', template=f) as g:
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
        write_singleton_reads(singleton_matrix, '17', g)

    # print(seqs)
    # consensus_seq=get_consensus_seq2(position_matrix)
    # print(position_matrix['CATGGCGAGCAT'])
    # print(umis['CATGGCGAGCCT'].count)
    # print(umis['CATGGCGAGCCT'].centroid)
    # for k in consensus_seq:
    #     if not consensus_seq[k]:
    #         print(None)
    #     else:
    #         print('@'+k+'_'+str(umis[k].count)+'_'+str(consensus_seq[k].start_pos)+'\n'+consensus_seq[k].seq+'\n+\n'+consensus_seq[k].qual)
    # fasta.close()
    # for k in singleton_matrix:
    #     print(k, umis[k].count, singleton_matrix[k].alignment.qname)
    # print(position_matrix['GCACCCGCGCAC'])
    # print(position_matrix['GCACCCGGGCCC'])


if __name__ == '__main__':
    main(sys.argv[1])
