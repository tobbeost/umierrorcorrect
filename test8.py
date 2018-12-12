#!/usr/bin/env python3
from group_new import readBam,read_bam_from_bed
from testgroup5 import umi_cluster,cluster_barcodes,get_connected_components,merge_clusters
from get_cons_dict2 import get_cons_dict, get_all_consensus,write_singleton_reads,consensus_read
import sys
import pysam
from umi_error_correct2 import cluster_umis_on_position
from collections import Counter
def umi_cluster(umi_dict):
    edit_distance_threshold=1
    adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
    clusters=get_connected_components(umi_dict,adj_matrix)
    umis=merge_clusters(umi_dict,clusters)
    return(umis)

def get_reads(bamfile,contig,start,end):
    with pysam.AlignmentFile(bamfile,'rb') as f:
        umis=Counter()
        reads=f.fetch(contig,start,end+1)
        for read in reads:
            barcode=read.qname.split(':')[-1]
            umis[barcode]+=1
            print(read.pos)
    return(umis)    

def main(bamfilename):
    regions,ends=cluster_umis_on_position(bamfilename)
    region=regions['2'][33141405]
    end=ends['2'][33141405]
    #umis=umi_cluster(region)
    print('2',33141405,end)
    rumis=get_reads(bamfilename,'2',33141405,end)
    #position_matrix,singleton_matrix=get_cons_dict(bamfilename,umis,'2',33141405,end)
    #consensus_seq=get_all_consensus(position_matrix,umis)
    #print(consensus_seq)
    #print(len(region))
    #print(len(umis))
    #print(len(rumis))
    #print(umis['ATGGGGGGGGGG'])
if __name__=='__main__':
    main(sys.argv[1])
