#!/usr/bin/env python3

from group import readBam
from testgroup5 import umi_cluster,cluster_barcodes,get_connected_components,merge_clusters
import sys
import pickle

def main(filename):
    regions=readBam(filename)
    print('cluster barcodes')
    edit_distance_threshold=1
    for chrx in regions:
        for r in regions[chrx]:
            adj_matrix=cluster_barcodes(regions[chrx][r],edit_distance_threshold)
            clusters=get_connected_components(regions[chrx][r],adj_matrix)
    #print(adj_matrix)
            umis=merge_clusters(regions[chrx][r],clusters)
            print(chrx,r,len(umis))
    #with open('test_new.pickle','wb') as f:
    #    pickle.dump(umis,f)


if __name__=='__main__':
    main(sys.argv[1])
