#!/usr/bin/env python
from group import readBam,read_bam_from_bed
from testgroup5 import umi_cluster,cluster_barcodes,get_connected_components,merge_clusters
from get_cons_dict2 import get_cons_dict, get_all_consensus,write_singleton_reads,consensus_read
import sys
import pysam
from multiprocessing import Pool

def umi_cluster_worker(args):
    umi_dict, region_id, edit_distance_threshold=args
    adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
    clusters=get_connected_components(umi_dict,adj_matrix)
    umis=merge_clusters(umi_dict,clusters)
    return(umis)

    '''Cluster UMIs based on UMI sequence for one region'''


def cluster_consensus_worker(args):
    umi_dict, region_id, edit_distance_threshold,bamfilename=args
    adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
    clusters=get_connected_components(umi_dict,adj_matrix)
    umis=merge_clusters(umi_dict,clusters)
    position_matrix,singleton_matrix=get_cons_dict(bamfilename,umis)
    consensus_seq=get_all_consensus(position_matrix,umis)
    outfilename='tmp/tmp_{}.bam'.format(region_id)
    #with pysam.AlignmentFile(bamfilename,'rb') as f, pysam.AlignmentFile(outfilename,'wb',template=f) as g:

    return(umis)
#def consensus_worker():

def cluster_umis_all_regions(regions,edit_distance_threshold,num_cpus):
    argvec=[]
    i=0
    for contig in sorted(regions):
        for pos in sorted(regions[contig]):
            argvec.append((regions[contig][pos],i,edit_distance_threshold))
            i+=1
        
    p=Pool(int(num_cpus))
    #print(argvec)
    p.map(umi_cluster_worker,argvec)

def cluster_umis_on_position(bamfilename,bedfilename=None):
    position_threshold=20
    group_method='fromBed'
    if group_method=='fromBed':
        regions=read_bam_from_bed(bamfilename,bedfilename,position_threshold)
    else:
        regions=readBam(bamfilename,position_threshold)
    for chrx in regions:
        regions2=regions[chrx]
        for rr in regions2:
            print(chrx,rr,regions2[rr].most_common(10)) 
    return(regions)

def main(bamfilename,bedfilename):
    regions=cluster_umis_on_position(bamfilename,bedfilename)
    for chrx in regions:
        for r in regions[chrx]:
            print(chrx,r,len(regions[chrx][r]),regions[chrx][r].most_common(1))
    edit_distance_threshold = 1
    num_cpus=40
    cluster_umis_all_regions(regions,edit_distance_threshold,num_cpus)
    ##print(regions)

if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])


