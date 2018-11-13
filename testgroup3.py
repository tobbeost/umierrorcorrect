#!/usr/bin/env python3
from collections import Counter
import time
import itertools
import sys
import pickle

class umi_cluster:
    def __init__(self,name,count):
        self.centroid=name
        self.count=count
    def add_count(self,newcount):
        self.count+=newcount
    def change_centroid(self,newname):
        self.centroid=newname

def hamming_distance(a, b):
    """Returns the Hamming distance between two strings of equal length"""
    assert len(a) == len(b) 
    return sum(i != j for i , j in zip(a, b))

def cluster_barcodes(barcodedict):
    adj_matrix={}
    comb=itertools.combinations(barcodedict.keys(),2)
    #print(comb)
    for a,b in comb:
        if hamming_distance(a,b) <= 1:
            if barcodedict[a] >= barcodedict[b]:
                if a not in adj_matrix:
                    adj_matrix[a]=[]
                adj_matrix[a].append(b)
            else:
                if b not in adj_matrix:
                    adj_matrix[b]=[]
                adj_matrix[b].append(a)
    return(adj_matrix)

def merge_clusters(barcodedict,adj_matrix):
    umis={}
    for name,count in barcodedict.items():
        umis[name]=umi_cluster(name,count)
    for centroid in adj_matrix:
        for neighbor in adj_matrix[centroid]:
            umis[centroid].add_count(umis[neighbor].count)
        for neighbor in adj_matrix[centroid]:
            umis[neighbor]=umis[centroid]
    #print(umis)
    return(umis)

def main():
    #print(hamming_distance('ATTAA','ACTAA'))
    #print(hamming_distance('ATAA','ACTAA'))
    with open('/home/xsteto/tmp/umierrorcorrect/test.pickle','rb') as f:
        regions=pickle.load(f)
        #adj_matrix=cluster_barcodes(regions['2'][33141588])
        adj_matrix=cluster_barcodes(regions['17'][7577495])
        print(adj_matrix)
        #umis=merge_clusters(regions['2'][33141588],adj_matrix)
        umis=merge_clusters(regions['17'][7577495],adj_matrix)
        print(umis)
        for umi in umis:
            print(umi,umis[umi].centroid,umis[umi].count)
    with open('/home/xsteto/umierrorcorrect/umi.pickle','wb') as g:
       pickle.dump(umis,g)
if __name__=='__main__':
    main()
