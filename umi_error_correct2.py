#!/usr/bin/env python
from group_new import readBam,read_bam_from_bed
from testgroup5 import umi_cluster,cluster_barcodes,get_connected_components,merge_clusters
from get_cons_dict2 import get_cons_dict, get_all_consensus,write_singleton_reads,consensus_read
import sys
import os
import pysam
from multiprocessing import Pool,cpu_count
import argparse

def parseArgs():
    parser=argparse.ArgumentParser(description="Pipeline for analyzing  barcoded amplicon sequencing data with Unique molecular identifiers (UMI)")
    parser.add_argument('-o', '--output_path',dest='output_path', help='Path to the output directory, required', required=True)
    parser.add_argument('-b','--bam', dest='bam_file',help='Path to BAM-file')
    parser.add_argument('-d','--edit_distance', dest='edit_distance_threshold', help="Edit distance threshold for UMI clustering, [default = %(default)s]",default=1)
    parser.add_argument('-p','--position_threshold',dest='position_threshold',help='Position threshold for grouping by position [default = %(default)s]',default=10)
    parser.add_argument('-singletons','--include_singletons', dest='include_singletons',action='store_true',help='Include this flag if singleton reads should be included in the output consensus read bam file. Note that the singletons will not be error corrected')
    parser.add_argument('-t','--num_threads',dest='num_threads',help='Number of threads to run the program on. If excluded, the number of cpus are automatically detected')
    args=parser.parse_args(sys.argv[1:])
    return(args)

#def umi_cluster_worker(args):
#    umi_dict, region_id, edit_distance_threshold=args
#    adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
#    clusters=get_connected_components(umi_dict,adj_matrix)
#    umis=merge_clusters(umi_dict,clusters)
#    return(umis)
#
#    '''Cluster UMIs based on UMI sequence for one region'''


def cluster_consensus_worker(args):
    umi_dict, tmpfilename, contig, start, end,  edit_distance_threshold,bamfilename,include_singletons=args
    adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
    clusters=get_connected_components(umi_dict,adj_matrix)
    umis=merge_clusters(umi_dict,clusters)
    position_matrix,singleton_matrix=get_cons_dict(bamfilename,umis,contig,start,end,include_singletons)
    consensus_seq=get_all_consensus(position_matrix,umis,contig)
    outfilename=tmpfilename
    with pysam.AlignmentFile(bamfilename,'rb') as f, pysam.AlignmentFile(outfilename,'wb',template=f) as g:
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
        if include_singletons:
            write_singleton_reads(singleton_matrix,contig,g)
    
#def consensus_worker():
#def merge_bams(output_path,bamfilelist):
#    print(output_path+'/consensus_reads.bam')
#    with open(output_path+'/consensus_reads.bam','wb') as g:
#        filename=bamfilelist[0]
#        print(filename)
#        with open(filename,'rb') as f:
#            for line in f:
#                g.write(line)
#        for filename in bamfilelist[1:]:
#            print(filename)
#            with open(filename,'rb') as f:
#                for line in f:
#                    if not line.startswith(b'@'): 
#                    #print(line)
#                        g.write(line)

def merge_bams(output_path,bamfilelist):
    with pysam.AlignmentFile(bamfilelist[0],'rb') as f,  pysam.AlignmentFile(output_path+'/consensus_reads.bam','wb',template=f) as g:
        for line in f:
            g.write(line)
        if len(bamfilelist)>1:
            for filename in bamfilelist[1:]:
                print(filename)
                with pysam.AlignmentFile(filename,'rb') as f1:
                    for line in f1:
                        g.write(line)
    for filename in bamfilelist:
        os.remove(filename)

def cluster_umis_all_regions(regions,ends,edit_distance_threshold,bamfilename,output_path,include_singletons,num_cpus):
    argvec=[]
    bamfilelist=[]
    i=0
    for contig in regions:
        for pos in regions[contig]:
            tmpfilename='{}/tmp_{}.bam'.format(output_path,i)
            argvec.append((regions[contig][pos],tmpfilename,contig,pos,ends[contig][pos],edit_distance_threshold,bamfilename,include_singletons))
            bamfilelist.append('{}/tmp_{}.bam'.format(output_path,i))
            i+=1
        
    p=Pool(int(num_cpus))
    #print(argvec)
    p.map(cluster_consensus_worker,argvec)
    return(bamfilelist)

def cluster_umis_on_position(bamfilename,position_threshold):
    position_threshold=20
    group_method='fromBed'
    group_method='automatic'
    if group_method=='fromBed':
        regions=read_bam_from_bed(bamfilename,bedfilename,position_threshold)
    else:
        regions,ends=readBam(bamfilename,position_threshold)
    #for chrx in regions:
    #    regions2=regions[chrx]
    #    for rr in regions2:
    #        print(chrx,rr,regions2[rr].most_common(10)) 
    return(regions,ends)

def main(args):
    regions,ends=cluster_umis_on_position(args.bam_file,args.position_threshold)
    for chrx in regions:
        for r in regions[chrx]:
            print(chrx,r,ends[chrx][r],len(regions[chrx][r]),regions[chrx][r].most_common(1))
    edit_distance_threshold = args.edit_distance_threshold
    if args.num_threads:
        num_cpus=args.num_threads
    else:
        num_cpus=cpu_count()
    print(num_cpus)
    bamfilelist=cluster_umis_all_regions(regions,ends,edit_distance_threshold,args.bam_file,args.output_path,args.include_singletons,num_cpus)
    merge_bams(args.output_path,bamfilelist)
    ##print(regions)

if __name__=='__main__':
    args=parseArgs()
    main(args)


