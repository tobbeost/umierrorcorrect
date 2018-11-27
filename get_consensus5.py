#!/usr/bin/env python3

import sys
import pickle
import pysam
from testgroup3 import umi_cluster
from collections import Counter
from math import log10

def get_reference_sequence(fasta,chrx,start,stop):
    ref=fasta.fetch(chrx,start,stop).upper()
    return(ref)
def get_consensus_bam(position_matrix,fasta,freq_threshold):
    consensus=position_matrix['GCACCCGCGCCC']
    consensusseq={}
    for position in consensus:
        if len(consensus[position])==1: #only one option
            consensusseq[position]=list(consensus[position].keys())[0]
        else:
            cons_allele = max(consensus[position], key = consensus[position].get) #get the allele with highest count
            cons_percent = (consensus[position][cons_allele]/sum(consensus[position].values())) * 100
            if cons_percent>=freq_threshold:
                consensusseq[position]=cons_allele
            else:
                consensusseq=None #failed to generate consensus read
                break
    return(consensusseq)

def get_consensus_seq2(position_matrix):
    consensus_sequences={}
    for barcode in position_matrix:
        seq=''
        qual=''
        for pos in position_matrix[barcode]:
            cons_base,cons_qual=calc_consensus_probabilities(position_matrix[barcode][pos])
            seq=seq+cons_base
            qual=qual+get_ascii(cons_qual)
        consensus_sequences[barcode]=seq+' '+qual
    return(consensus_sequences)
#def get_cigar_string(cons,ref):
#    cons_sequence=''.join(cons.values())
#    if 'D' in cons_sequence or 'I' in cons_sequence:
#
#    else:
#        if cons_sequence==ref:
#            return('{}M'.format(len(cons)))
#        else:
#        q
def get_phred(character):
    '''Returns the numeric value of ASCII character associated with phred score (offset 33)'''
    value=ord(character)-33
    return(value)

def get_ascii(value):
    ascii_letter=chr(value+33)
    return(ascii_letter)

def calc_consensus(base,cons_pos):
    prod=1
    for nucl in cons_pos:
        if nucl in base:
            for phred in cons_pos[nucl]:
                prod=prod*(1-(10**(-phred/10)))
        if nucl not in base and nucl in 'ATCG':
            for phred in cons_pos[nucl]:
                prod=prod*(10**(-phred/10))
    return(prod)


def calc_consensus_probabilities(cons_pos):
    p={base:calc_consensus(base,cons_pos) for base in 'ATCG'}
    denom=sum(p.values())
    probs={base:p[base]/denom for base in 'ATCG'}
    cons_base=max(probs, key=probs.get)
    if probs[cons_base]==1:
        cons_phred=60
    else:
        cons_phred=round(-10*log10(1-probs[cons_base]))
        if cons_phred > 60:
            cons_phred=60
    return(cons_base,cons_phred)

def get_consensus(bamfilename,umis,family_sizes,ref_seq):
    consensus_sequence={}
    for fs in family_sizes:
        consensus_sequence[fs]={}
    position_matrix={}
    singleton_matrix={}
    with pysam.AlignmentFile(bamfilename,'rb') as f:
        alignment=f.pileup('17',7577497,7577600,max_depth=1000000)
        n=0
        clustlist=['GCACCCGCGCCC','GCACCCGCGCAC','GCACCCGGGCCC']
        for pileupcolumn in alignment:
            pos=pileupcolumn.pos
            #ref_base=ref_seq[pos-7577497]
            #print(pos)
            for read in pileupcolumn.pileups:
                barcode=read.alignment.qname.split(':')[-1]
                if barcode in clustlist and pos==7577513:
                    n+=1
                cluster=umis[barcode].centroid
                cluster_size=umis[barcode].count
                if cluster_size > 1:    
                    if not read.is_del and not read.indel:
                        q=get_phred(read.alignment.qual[read.query_position])
                        al=read.alignment.query_sequence[read.query_position]
                        #if cluster in 'GCACCCGCGCCC' and pos==7577497:
                        #    print(allele)
                        #    n+=1
                    elif read.indel>0: #insertion
                        #print(read.indel)
                        al='I'
                        q=get_phred(read.alignment.qual[read.query_position])
                        #al=read.alignment.query_sequence[read.query_position:read.query_position+abs(read.indel)+1]
                    elif read.indel<0: #deletion
                        #print(read.indel)
                        #al=read.alignment.query_sequence[read.query_position+1]
                    #ref_base=ref_seq[(pos-7577497):(pos-7577497)+abs(read.indel)]
                        al='D'
                        q=get_phred(read.alignment.qual[read.query_position])
                    if cluster not in position_matrix:
                        position_matrix[cluster]={}
                    if pos not in position_matrix[cluster]:
                        position_matrix[cluster][pos]={}
                    #allele=(ref_base,al)
                    allele=al
                    if allele not in position_matrix[cluster][pos]:
                        position_matrix[cluster][pos][allele]=[]
                    position_matrix[cluster][pos][allele].append(q)
                else:
                    if cluster not in singleton_matrix:
                        singleton_matrix[cluster]=read


    print(n)
    return(position_matrix,singleton_matrix)



def main(bamfilename):
    with open('/home/xsteto/tmp/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    family_sizes=[0,1,2,3,4,5,7,10,20,30]
    fasta=pysam.FastaFile('/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta')
    ref_seq=get_reference_sequence(fasta,'17',7577497,7577800)

    position_matrix,singleton_matrix=get_consensus(bamfilename,umis,family_sizes,ref_seq)
    consensus_seq=get_consensus_seq2(position_matrix)
    print(consensus_seq['GCACCCGCGCCC'])
    for k in consensus_seq:
        print(k,umis[k].count,consensus_seq[k])
    fasta.close()
    #for k in singleton_matrix:
    #    print(k, umis[k].count, singleton_matrix[k].alignment.qname)

    #print(position_matrix['GCACCCGCGCAC'])
    #print(position_matrix['GCACCCGGGCCC'])
if __name__=='__main__':
    main(sys.argv[1])
