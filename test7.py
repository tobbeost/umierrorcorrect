#!/usr/bin/env python3
import pysam
from get_cons_dict import *
from testgroup5 import *
import pickle

def getConsensus4(group_seqs,indel_freq_threshold,umi_info):
    '''Takes a list of pysam entries (rows in the BAM file) as input and generates a consensus sequence.'''
    consensus={}
    is_indel=False
    template_read=group_seqs[0]
    for read in group_seqs:
        if 'I' not in read.cigarstring and 'D' not in read.cigarstring: #no indel in read
            sequence=read.seq
            qual=read.qual
            for qpos,refpos in read.get_aligned_pairs(matches_only=True):
                base=sequence[qpos]
                if refpos not in consensus:
                    consensus[refpos]={}
                if base not in consensus[refpos]:
                    consensus[refpos][base]=[]
                consensus[refpos][base].append(get_phred(qual[qpos]))
        else: #indel at next position
            positions=read.get_aligned_pairs(matches_only=True)
            q,ref=positions[0]
            for qpos,refpos in positions[1:]:
                if not qpos==q+1:
                    print('insertion',(q+1),qpos,read.seq[q+1:qpos])
                    allele=read.seq[q+1:qpos]
                    if refpos not in consensus:
                        consensus[refpos]={}
                    if 'I' not in consensus[refpos]:
                        consensus[refpos]['I']={}
                    if allele not in consensus[refpos]['I']:
                        consensus[refpos]['I'][allele]=0
                    consensus[refpos]['I'][allele]+=1
                elif not refpos==ref+1:
                    print('deletion',ref+1,refpos)
                    dellength=refpos-(ref+1)
                    if refpos not in consensus:
                        consensus[refpos]={}
                    if 'D' not in consensus[refpos]:
                        consensus[refpos]['D']={}
                    if dellength not in consensus[refpos]['D']:
                        consensus[refpos]['D'][dellength]=0
                    consensus[refpos]['D'][dellength]+=1
                else:
                    sequence=read.seq
                    qual=read.qual
                    base=sequence[qpos]
                    if refpos not in consensus:
                        consensus[refpos]={}
                    if base not in consensus[refpos]:
                        consensus[refpos][base]=[]
                    consensus[refpos][base].append(get_phred(qual[qpos]))
                q=qpos
                ref=refpos
    return(consensus)

def get_read(bc):
    with open('/home/xsteto/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    bamfilename='output.sorted.bam'
    position_matrix,singleton_matrix=get_cons_dict(bamfilename,umis)
    read=getConsensus3(position_matrix[bc],50.0,umis[bc])
    return(read)
def main():
    bc='CGTAGTTGTCCT'
    read=get_read(bc)
    print(read.get_cigar())
    print(read.cigarstring)
if __name__=='__main__':
    main()
