#!/usr/bin/env python3

from __future__ import division
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
        a.tags = (("NM", self.nmtag), ("RG", "L1"))
        f.write(a)


def get_reference_sequence(fasta, chrx, start, stop):
    '''Returns the fasta sequence of the reference for a given region'''
    ref = fasta.fetch(chrx, start, stop).upper()
    return(ref)


def get_most_common_allele(cons_pos):
    '''Calculate the allele frequencies at one position and returns the 
       allele with the highest frequency.'''
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


def get_phred(character):
    '''Get the numeric value of ASCII character associated with phred score (offset 33)'''
    value = ord(character) - 33
    return(value)


def get_ascii(value):
    '''Get the ascii character for a given phred score (offset 33)'''
    ascii_letter = chr(value + 33)
    return(ascii_letter)


def calc_consensus(base, cons_pos):
    '''Function for calculating the combined score for a base at a position'''
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
    '''Function for calculating the probability of consensus at a given position
       Return the base with the highest probability'''
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


def getConsensus3(group_seqs, contig, regionid, indel_freq_threshold, umi_info, consensus_freq_threshold):
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
                        dellength = refpos - (ref+1) #get the length of the deletion
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
    
    if len(consensus) > 0:
        #generate the consensus sequence
        consensus_sorted = sorted(consensus)
        consread = consensus_read(contig, regionid, consensus_sorted[0], umi_info.centroid, umi_info.count)
        add_consensus = True
        skippos = [] #if position is del
        prevpos = consensus_sorted[0] - 1
        for pos in sorted(consensus_sorted):
            if pos not in skippos:
                if not pos == prevpos + 1:
                    for i in range(prevpos+1,pos):
                        if i not in skippos:
                            consread.add_base('N', get_ascii(0))
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
                            if dellength > 1:
                                for i in range(1,dellength):
                                    skippos.append(pos + i)
                        else:
                            consread.add_base('N', get_ascii(0))
                            add_consensus = False
                    elif percent >= indel_freq_threshold:
                        cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                        consread.add_base(cons_base, get_ascii(cons_qual))
                    else:
                        consread.add_base('N', get_ascii(0))
                        add_consensus = False
                else:
                    #no indel
                    cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                    if consensus_freq_threshold: #test if not None
                        if len(consensus[pos]) == 1:  #100%
                            consread.add_base(cons_base, get_ascii(cons_qual))
                        else:
                            if cons_base not in consensus[pos]:
                                print(cons_base+" not in consensus[pos] "+str(pos), consensus[pos])
                            else:
                                percent = (len(consensus[pos][cons_base]) / len(group_seqs))*100.0
                            if percent >= consensus_freq_threshold: #consensus frequency above threshold
                                consread.add_base(cons_base, get_ascii(cons_qual))
                            else:
                                consread.add_base('N', get_ascii(0))
                                add_consensus = False
                    else:
                        consread.add_base(cons_base, get_ascii(cons_qual))
            prevpos=pos

        if add_consensus:
            return(consread)
        else:
            return(None)
    else:
        return(None)


def get_all_consensus(position_matrix, umis, contig, regionid, indel_frequency_cutoff, consensus_frequency_cutoff):
    '''Get the consensus sequences for all umis'''
    consensuses = {}
    for umi in position_matrix:
        consensuses[umi] = getConsensus3(position_matrix[umi], contig, regionid, 
                                         indel_frequency_cutoff, umis[umi],
                                         consensus_frequency_cutoff)
    return(consensuses)


def get_cons_dict(bamfilename, umis, contig, start, end, include_singletons):
    position_matrix = {}
    singleton_matrix = {}
    with pysam.AlignmentFile(bamfilename, 'rb') as f:
        alignment = f.fetch(contig, start, end)
        for read in alignment:
            barcode = read.qname.split(':')[-1]
            pos = read.pos
            if pos >= start and pos <= end:
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


def write_singleton_reads(singleton_matrix, region_id, g):
    for umi, read in singleton_matrix.items():
        read.query_name = 'Singleton_read_{}_{}_Count=1'.format(region_id, umi)
        g.write(read)


def main(bamfilename):
    with open('/home/xsteto/umierrorcorrect/umi.pickle', 'rb') as f:
        umis = pickle.load(f)
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, contig, start, end)
    consensus_seq = get_all_consensus(position_matrix, umis, contig)
    with pysam.AlignmentFile(bamfilename, 'rb') as f, pysam.AlignmentFile('consensus_out.bam', 'wb', template=f) as g:
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
        write_singleton_reads(singleton_matrix, '17', g)



if __name__ == '__main__':
    main(sys.argv[1])
