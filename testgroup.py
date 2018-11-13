#!/usr/bin/env python3
from collections import Counter
import time
import itertools
import sys
import pickle
def call_counter(func):
    def helper(*args, **kwargs):
        helper.calls += 1
        return func(*args, **kwargs)
    helper.calls = 0
    helper.__name__= func.__name__
    return helper
def memoize(func):
    mem = {}
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in mem:
            mem[key] = func(*args, **kwargs)
        return mem[key]
    return memoizer
@call_counter
@memoize
def levenshtein(s, t):
    if s == "":
        return len(t)
    if t == "":
        return len(s)
    if s[-1] == t[-1]:
        cost = 0
    else:
        cost = 1

    res = min([levenshtein(s[:-1], t)+1,
               levenshtein(s, t[:-1])+1,
               levenshtein(s[:-1], t[:-1]) + cost])
    return res



def cluster_barcodes(barcodedict):
    adj_matrix={}
    comb=itertools.combinations(barcodedict.keys(),2)
    #print(comb)
    for a,b in comb:
        if levenshtein(a,b) <= 1:
            if barcodedict[a] >= barcodedict[b]:
                if a not in adj_matrix:
                    adj_matrix[a]=[]
                adj_matrix[a].append(b)
            else:
                if b not in adj_matrix:
                    adj_matrix[b]=[]
                adj_matrix[b].append(a)
    return(adj_matrix)

def main():
    with open('test.pickle','rb') as f:
        regions=pickle.load(f)
        adj_matrix=cluster_barcodes(regions['17'][7577495])
        print(adj_matrix)
if __name__=='__main__':
    main()
