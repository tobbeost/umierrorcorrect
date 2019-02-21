#/usr/bin/env python3
import pickle
from src.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from src.get_consensus import get_cons_dict, getConsensus3, write_singleton_reads, get_reference_sequence,get_phred
import time

def cluster():
    with open('/home/xsteto/tmp/umierrorcorrect/test.pickle','rb') as f:
        regions=pickle.load(f)
    umi_dict=regions['17'][7577495]
    edit_distance_threshold=1
    adj_matrix = cluster_barcodes(umi_dict, edit_distance_threshold)
    clusters=get_connected_components(umi_dict,adj_matrix)
    umis=merge_clusters(umi_dict,clusters)
    li=['CATGGCGAGCAT', 'CATGGCGAGCCT', 'CATGGCTAGCAT']
    for u in umis:
        if u in li:
            print(u,umis[u].centroid,umis[u].count)
    return(umis)

def consensus(umis,bamfilename):
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, '17', 7577495, 7577600, True)
    consensuses = {}
    for umi in position_matrix:
        consensuses[umi] = getConsensus3(position_matrix[umi], 17, 1, 75.0, umis[umi])
    
    return(consensuses)

def getConsensus2(group_seqs,contig,regionid,indel_freq_threshold,umi_info):
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
    return(consensus)

def consensus_test(umis,bamfilename):
    start=time.time()
    u2={'TTTCACAACGCA': umis['TTTCACAACGCA']}
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, '17', 7577495, 7577600, True)
    #seqs=position_matrix['TTTCACAACGCA']
    #seqs=position_matrix
    #umis_sel=umis['TTTCACAACGCA']
    #umis_sel=umis
    consensuses={}
    for umi in position_matrix:
        seqs=position_matrix[umi]
        umis_sel=umis[umi]
        consensuses[umi] = getConsensus2(seqs, 17, 1, 75.0, umis_sel)
    end=time.time()
    print(end-start)
    return(consensuses)

def main():
    umis=cluster()
    bamfilename='/home/xsteto/testuec/output.sorted.bam'
    cons=consensus_test(umis,bamfilename)
    #print(cons)

if __name__=='__main__':
    main()

