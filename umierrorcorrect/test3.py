#!/usr/bin/env python3
import pickle
from umierrorcorrect.src.umi_cluster import cluster_barcodes, get_connected_components, merge_clusters
from umierrorcorrect.src.get_consensus import get_cons_dict, getConsensus3, write_singleton_reads, get_reference_sequence,get_phred, consensus_read, calc_consensus_probabilities,get_ascii,get_most_common_allele
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


def generate_consensus(consensus,group_seqs,contig,regionid,indel_freq_threshold,umi_info,consensus_freq_threshold):
    if len(consensus) > 0:
        consensus_sorted = sorted(consensus)
        consread = consensus_read(contig, regionid, consensus_sorted[0], umi_info.centroid, umi_info.count)
        add_consensus = True
        for pos in sorted(consensus_sorted):
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
                    else:
                        consread.add_base('N', get_ascii(0))
                        add_consensus = False
                elif percent >= indel_freq_threshold:
                    cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                    consread.add_base(cons_base, get_ascii(cons_qual))
                else:
                    consread.add_base('N', get_ascii(0))
            else:
                #no indel
                cons_base, cons_qual = calc_consensus_probabilities(consensus[pos])
                if consensus_freq_threshold:
                    if len(consensus[pos]) == 1:  #100%
                        consread.add_base(cons_base, get_ascii(cons_qual))
                    else:
                        percent = (len(consensus[pos][cons_base]) / len(group_seqs))*100.0
                        if percent >= consensus_freq_threshold:
                            consread.add_base(cons_base, get_ascii(cons_qual))
                        else:
                            consread.add_base('N', get_ascii(0))
                            add_consensus = False
                else:
                    consread.add_base(cons_base, get_ascii(cons_qual))
        if add_consensus:
            return(consread)
        else:
            print(umi_info.centroid,umi_info.count)
            return(None)
    else:
        return(None)

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
    #consreads={}
    #for umi in consensus:
    #    seqs=position_matrix[umi]

    end=time.time()
    print(end-start)
    return(consensuses)

def main():
    umis=cluster()
    bamfilename='/home/xsteto/testuec/output.sorted.bam'
    consensus_freq_threshold=60.0
    #consensus_freq_threshold=None
    #cons=consensus_test(umis,bamfilename)
    #with open('consensus_matrix.pickle','wb') as f:
    #    pickle.dump(cons,f)
    position_matrix, singleton_matrix = get_cons_dict(bamfilename, umis, '17', 7577495, 7577600, True)
    with open('consensus_matrix.pickle','rb') as f:
        cons=pickle.load(f)
    start=time.time()
    consensuses={}
    for umi in cons:
        seqs=cons[umi]
        umis_sel=umis[umi]
        group_seqs=position_matrix[umi]
        consensuses[umi] = generate_consensus(seqs,group_seqs, 17, 1, 75.0, umis_sel,consensus_freq_threshold)
    end=time.time()
    print(end-start)
    print(consensuses['CATGGCGAGCAT'].seq)

if __name__=='__main__':
    main()

