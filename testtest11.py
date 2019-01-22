#!/usr/bin/env python3
import sys
from test13 import get_cons_info,write_consensus
from umi_error_correct4 import cluster_umis_on_position,cluster_consensus_worker
from get_regions_from_bed2 import read_bed,sort_regions,merge_regions,get_annotation
from get_cons_dict2 import get_cons_dict, get_all_consensus,write_singleton_reads,consensus_read,get_reference_sequence
from testgroup5 import umi_cluster,cluster_barcodes,get_connected_components,merge_clusters
import pysam
import pickle

def cluster_consensus_worker(args):
    umi_dict, tmpfilename, regionid, contig, start, end,  edit_distance_threshold,bamfilename,include_singletons,annotations,fasta=args
    adj_matrix=cluster_barcodes(umi_dict,edit_distance_threshold)
    clusters=get_connected_components(umi_dict,adj_matrix)
    umis=merge_clusters(umi_dict,clusters)
    position_matrix,singleton_matrix=get_cons_dict(bamfilename,umis,contig,start,end,True) #include_singletons=True
    consensus_seq=get_all_consensus(position_matrix,umis,contig,regionid)
    outfilename=tmpfilename
    with pysam.AlignmentFile(bamfilename,'rb') as f, pysam.AlignmentFile(outfilename,'wb',template=f) as g:
        for cons_read in consensus_seq.values():
            if cons_read:
                cons_read.write_to_bam(g)
        if include_singletons:
            write_singleton_reads(singleton_matrix,contig,g)
    with open('testnew5/cons.pickle','wb') as g:
        pickle.dump(consensus_seq,g)
    #with open('testnew5/single.pickle','wb') as g:
    #    pickle.dump(singleton_matrix,g)
    cons=get_cons_info(consensus_seq,singleton_matrix)

    consfilename=outfilename.rstrip('.bam')+'.cons'
    statfilename=outfilename.rstrip('.bam')+'.hist'
    startpos=min(list(cons.keys()))
    endpos=max(list(cons.keys()))+1
    with pysam.FastaFile(fasta) as f:
        ref_seq=get_reference_sequence(f,contig,startpos,endpos)
    with open(consfilename,'w') as g:
        write_consensus(g,cons,ref_seq,startpos,contig,annotations['3'],False)
    with open(statfilename,'w') as g2:
        regionname='{}:{}-{}'.format(contig,start,end)
        g2.write('\t'.join([str(regionid),regionname,'singletons: '+str(len(singleton_matrix))])+'\n')


def main():
    bam_file='output.sorted.bam'
    regions,ends=cluster_umis_on_position(bam_file,10)
    for contig in regions:
        print(contig,regions[contig].keys())
    edit_distance_threshold=1
    contig='3'
    pos=178952035
    tmpfilename='testnew7/tmp_9.bam'
    bedfile='/home/xsteto/Cons_depth_estimation/assay_regions.txt'
    bedregions=read_bed(bedfile)
    bedregions=sort_regions(bedregions)
    bedregions=merge_regions(bedregions,0)
    print(bedregions)
    fasta='/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta'
    cluster_consensus_worker((regions[contig][pos],tmpfilename,9,contig,pos,ends[contig][pos],edit_distance_threshold,bam_file,True,bedregions,fasta))

if __name__=='__main__':
    main()
