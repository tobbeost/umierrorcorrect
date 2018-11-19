#!/usr/bin/env python3

from group import readBam
from testgroup3 import umi_cluster,cluster_barcodes,merge_clusters
import sys
import pickle

def main(filename):
    regions=readBam(filename)
    print('cluster barcodes')
    adj_matrix=cluster_barcodes(regions['17'][7577495])
    print(adj_matrix)
    umis=merge_clusters(regions['17'][7577495],adj_matrix)
    with open('test_new.pickle','wb') as f:
        pickle.dump(umis,f)


if __name__=='__main__':
    main(sys.argv[1])
