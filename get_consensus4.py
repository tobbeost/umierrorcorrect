#!/usr/bin/env python3

import sys
import pickle
import pysam
from testgroup3 import umi_cluster
from collections import Counter

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

#def get_cigar_string(cons,ref):
#    cons_sequence=''.join(cons.values())
#    if 'D' in cons_sequence or 'I' in cons_sequence:
#
#    else:
#        if cons_sequence==ref:
#            return('{}M'.format(len(cons)))
#        else:
#        q
        
def get_consensus(bamfilename,umis,family_sizes,ref_seq):
    consensus_sequence={}
    for fs in family_sizes:
        consensus_sequence[fs]={}
    position_matrix={}
    with pysam.AlignmentFile(bamfilename,'rb') as f:
        alignment=f.pileup('17',7577497,7577600,max_depth=1000000)
        n=0
        clustlist=['GCACCCGCGCCC','GCACCCGCGCAC','GCACCCGGGCCC']
        for pileupcolumn in alignment:
            pos=pileupcolumn.pos
            ref_base=ref_seq[pos-7577497]
            print(pos)
            for read in pileupcolumn.pileups:
                barcode=read.alignment.qname.split(':')[-1]
                if barcode in clustlist and pos==7577513:
                    n+=1
                cluster=umis[barcode].centroid
                cluster_size=umis[barcode].count
                    
                if not read.is_del and not read.indel:
                    al=read.alignment.query_sequence[read.query_position]
                    #if cluster in 'GCACCCGCGCCC' and pos==7577497:
                    #    print(allele)
                    #    n+=1
                elif read.indel>0: #insertion
                    #print(read.indel)
                    al=read.alignment.query_sequence[read.query_position:read.query_position+abs(read.indel)+1]
                elif read.indel<0: #deletion
                    print(read.indel)
                    al=read.alignment.query_sequence[read.query_position+1]
                    ref_base=ref_seq[(pos-7577497):(pos-7577497)+abs(read.indel)]
                if cluster not in position_matrix:
                    position_matrix[cluster]={}
                if pos not in position_matrix[cluster]:
                    position_matrix[cluster][pos]={}
                allele=(ref_base,al)
                if allele not in position_matrix[cluster][pos]:
                    position_matrix[cluster][pos][allele]=0
                position_matrix[cluster][pos][allele]+=1
    print(n)
    return(position_matrix)



def main(bamfilename):
    with open('/home/xsteto/tmp/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    family_sizes=[0,1,2,3,4,5,7,10,20,30]
    fasta=pysam.FastaFile('/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta')
    ref_seq=get_reference_sequence(fasta,'17',7577497,7577800)

    position_matrix=get_consensus(bamfilename,umis,family_sizes,ref_seq)
    print(position_matrix['GCACCCGCGCCC'])
    print(len(position_matrix))
    fasta.close()
    #print(position_matrix['GCACCCGCGCAC'])
    #print(position_matrix['GCACCCGGGCCC'])
if __name__=='__main__':
    main(sys.argv[1])
