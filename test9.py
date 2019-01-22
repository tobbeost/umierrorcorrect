#!/usr/bin/env python3

import sys
import pickle
import pysam
from testgroup3 import umi_cluster
from collections import Counter
from math import log10
from itertools import groupby
from get_cons_dict2 import *


def write_to_cons(consensus_seq,singletons,filename,fsizes=[0,1,2,3,4,5,7,10,20,30]):
    cons={}
    for consensus_read in consensus_seq.values():
        pos=consensus_read.start_pos
        count=consensus_read.count
        if consensus_read.indel_read==0:
            for base in consensus_read.seq:
                if pos not in cons:
                    cons[pos]={}
                for fsize in fsizes:
                    if count>=fsize:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize][base]+=1
                pos+=1
        else:
            i=0
            cigar=consensus_read.cigarstring
            for base in consensus_read.seq:
                c=cigar[i]
                if pos not in cons:
                    cons[pos]={}
                if c=='0':
                    for fsize in fsizes:
                        if count>=fsize:
                            if fsize not in cons[pos]:
                                cons[pos][fsize]=Counter()
                            cons[pos][fsize][base]+=1
                    pos+=1
                    i+=1 
                elif c=='1':
                    for fsize in fsizes:
                        if count>=fsize:
                            if fsize not in cons[pos]:
                                cons[pos][fsize]=Counter()
                            cons[pos][fsize]['I']+=1
                    i+=1
                elif c=='2':
                    for fsize in fsizes:
                        if count>=fsize:
                            if fsize not in cons[pos]:
                                cons[pos][fsize]=Counter()
                            cons[pos][fsize]['D']+=1
                    pos+=1
                    i+=1
    for read in singletons.values():
        pos=read.pos
        if 'I' not in read.cigarstring and 'D' not in read.cigarstring:
            for base in read.seq:
                if pos not in cons:
                    cons[pos]={}
                for fsize in [0,1]:
                    if fsize not in cons[pos]:
                        cons[pos][fsize]=Counter()
                    cons[pos][fsize][base]+=1
            pos+=1
        else:
            cigar=''.join([str(a)*b for (a,b) in read.cigar])
            i=0
            for base in read.seq:
                c=cigar[i]
                if pos not in cons:
                    cons[pos]={}
                if c=='0':
                    for fsize in [0,1]:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize][base]+=1
                    pos+=1
                    i+=1
                elif c=='1':
                    for fsize in [0,1]:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize]['I']+=1
                    i+=1
                elif c=='2':
                    for fsize in [0,1]:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize]['D']+=1
                    pos+=1
                    i+=1 
                else:
                    i+=1
    print(cons)

                


def get_all_consensus(position_matrix,umis,contig):
    consensuses={}
    for umi in position_matrix:
        consensuses[umi]=getConsensus3(position_matrix[umi],contig,50.0,umis[umi])
    return(consensuses)
    

def get_cons_dict(bamfilename,umis,contig,start,end,include_singletons):
    #print('{}:{}-{}'.format(contig,start,end))
    position_matrix={}
    singleton_matrix={}
    with pysam.AlignmentFile(bamfilename,'rb') as f:
        alignment=f.fetch(contig,start,end)
        for read in alignment:
            barcode=read.qname.split(':')[-1]
            pos=read.pos
            if pos>=start and pos <=end:
            #if barcode not in ['ATGGGGGGGGGG','TCGGGGGGGGGG','TAGGGGGTGGGG','TAGGGGGGGGGG']:
            #    print('{}:{}-{}'.format(contig,start,end))
            #    for u in umis:
            #        print(u, umis[u].centroid, umis[u].count)
                cluster=umis[barcode].centroid
                cluster_size=umis[barcode].count
                if cluster_size>1:
                    if cluster not in position_matrix:
                        position_matrix[cluster]=[]
                    position_matrix[cluster].append(read)
                elif include_singletons:
                    if cluster not in singleton_matrix:
                        singleton_matrix[cluster]=read
    return(position_matrix,singleton_matrix)

def write_singleton_reads(singleton_matrix,contig,g):
    for umi,read in singleton_matrix.items():
        read.query_name='Singleton_read_{}_{}_Count=1'.format(contig,umi)
        g.write(read)

def main(bamfilename):
    with open('/home/xsteto/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    #family_sizes=[0,1,2,3,4,5,7,10,20,30]
    #fasta=pysam.FastaFile('/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta')
    #ref_seq=get_reference_sequence(fasta,'17',7577497,7577800)
    contig='17'
    start=7577495
    end=7577800
    position_matrix,singleton_matrix=get_cons_dict(bamfilename,umis,contig,start,end,True)
    consensus_seq=get_all_consensus(position_matrix,umis,contig)
    write_to_cons(consensus_seq,singleton_matrix,'out/consensus.cons')
    #with pysam.AlignmentFile(bamfilename,'rb') as f, pysam.AlignmentFile('consensus_out.bam','wb',template=f) as g:
    
    #    for cons_read in consensus_seq.values():
    #        if cons_read:
    #            cons_read.write_to_bam(g)
    #    write_singleton_reads(singleton_matrix,'17',g)

    #print(seqs)
    #consensus_seq=get_consensus_seq2(position_matrix)
    #print(position_matrix['CATGGCGAGCAT'])
    #print(umis['CATGGCGAGCCT'].count)
    #print(umis['CATGGCGAGCCT'].centroid)
    #for k in consensus_seq:
    #    if not consensus_seq[k]:
    #        print(None)
    #    else:
    #        print('@'+k+'_'+str(umis[k].count)+'_'+str(consensus_seq[k].start_pos)+'\n'+consensus_seq[k].seq+'\n+\n'+consensus_seq[k].qual)
    #fasta.close()
    #for k in singleton_matrix:
    #    print(k, umis[k].count, singleton_matrix[k].alignment.qname)

    #print(position_matrix['GCACCCGCGCAC'])
    #print(position_matrix['GCACCCGGGCCC'])

if __name__=='__main__':
    main(sys.argv[1])
