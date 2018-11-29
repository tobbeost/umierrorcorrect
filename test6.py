#!/usr/bin/env python3
from testgroup5 import *
def main():
    with open('/home/xsteto/umierrorcorrect/umi.pickle','rb') as f:
        umis=pickle.load(f)
    li=['CATGGCGAGCAT', 'CATGGCGAGCCT', 'CATGGCTAGCAT']
    for u in umis:
        if u in li:
            print(u,umis[u].centroid,umis[u].count)
    print('hej')
    for u in umis:
        if umis[u].centroid in li:
            print(u,umis[u].centroid,umis[u].count)

if __name__=='__main__':
    main()

