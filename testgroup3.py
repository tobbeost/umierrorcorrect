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

def create_substring_matrix(barcodedict,edit_distance_threshold):
    """Divide each barcode in two or three substrings of (approximately) equal length"""
    umi_length=len(list(barcodedict.keys())[0])
    if edit_distance_threshold <= 1:
        print('hello')
        s=round(umi_length/2)
        substr_dict1={}
        substr_dict2={}
        for barcode in barcodedict:
            sub1=barcode[:s]
            sub2=barcode[s:]
            if sub1 not in substr_dict1:
                substr_dict1[sub1]=[]
            if sub2 not in substr_dict2:
                substr_dict2[sub2]=[]
            substr_dict1[sub1].append(barcode)
            substr_dict2[sub2].append(barcode)
        return([substr_dict1,substr_dict2])
    if edit_distance_threshold == 2:
        s=round(umi_length/3)
        substr_dict1={}
        substr_dict2={}
        substr_dict3=[]
        for barcode in barcodedict:
            sub1=barcode[:s]
            sub2=barcode[s:2*s]
            sub3=barcode[2*s:]
            if sub1 not in substr_dict1:
                substr_dict1[sub1]=[]
            if sub2 not in substr_dict2:
                substr_dict2[sub2]=[]
            if sub3 not in substr_dict3:
                substr_dict3[sub3]=[]
            substr_dict1[sub1].append(barcode)
            substr_dict2[sub2].append(barcode)
            substr_dict3[sub3].append(barcode)
        return([substr_dict1,substr_dict2,substr_dict3])

def get_adj_matrix_from_substring(barcodedict,substrdictlist):
    '''A generator that generates combinations to test for edit distance'''
    umi_length=len(list(barcodedict.keys())[0])
    if len(substrdictlist)==2:
        s=round(umi_length/2)
        for barcode in barcodedict:
            neighbors=set()
            sub1=barcode[:s]
            neighbors=neighbors.union(substrdictlist[0][sub1])
            sub2=barcode[s:]
            neighbors=neighbors.union(substrdictlist[1][sub2])
            neighbors.remove(barcode)
            for neighbor in neighbors:
                yield barcode, neighbor
    if len(substrdictlist)==3:
        substr_dict1,substr_dict2,substr_dict3=substrdictlist
        s=round(umi_length/3)
        for barcode in barcodedict:
            neighbors=set()
            sub1=barcode[:s]
            neighbors=neighbors.union(substr_dict1[sub1])
            sub2=barcode[s:2*s]
            neighbors=neighbors.union(substr_dict2[sub2])
            sub3=barcode[2*s:]
            neighbors=neighbors.union(substr_dict3[sub3])
            neighbors.remove(barcode)
            for neighbor in neighbors:
                yield barcode, neighbor
        
def cluster_barcodes(barcodedict,edit_distance_threshold):
    adj_matrix={}
    if len(barcodedict)>30: 
        #compare substrings for speedup
        substring_matrix=create_substring_matrix(barcodedict,edit_distance_threshold)
        comb=get_adj_matrix_from_substring(barcodedict,substring_matrix)
    else:
        comb=itertools.combinations(barcodedict.keys(),2)
    #print(comb)
    #comb=itertools.combinations(barcodedict.keys(),2)
    for a,b in comb:
        if hamming_distance(a,b) <= edit_distance_threshold:
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
    edit_distance_threshold=1
    with open('/home/xsteto/tmp/umierrorcorrect/test.pickle','rb') as f:
        regions=pickle.load(f)
        #adj_matrix=cluster_barcodes(regions['2'][33141588])
        adj_matrix=cluster_barcodes(regions['17'][7577495],edit_distance_threshold)
        #print(adj_matrix)
        #umis=merge_clusters(regions['2'][33141588],adj_matrix)
        umis=merge_clusters(regions['17'][7577495],adj_matrix)
        #print(umis)
        for umi in umis:
            print(umi,umis[umi].centroid,umis[umi].count)
    with open('/home/xsteto/umierrorcorrect/umi.pickle','wb') as g:
       pickle.dump(umis,g)
if __name__=='__main__':
    main()
