#!/usr/bin/env python3
import sys
import pickle
import itertools
def main(filename):
    with open('/home/xsteto/tmp/umierrorcorrect/test.pickle','rb') as f:
        regions=pickle.load(f)
    with open(filename) as f:
        barcodematrix=regions['17'][7577495]
        umis={}
        for u in barcodematrix:
            umis[u]=0
        for line in itertools.combinations(barcodematrix.keys(),2):
            print(line)
        #    line=line.rstrip()
        #    barcodes=line.split()
        #    for bc in barcodes:
        #        umis[bc]+=1
    #for umi,count in umis.items():
    #    if count>1:
    #        print(umi,count)

if __name__=='__main__':
    main(sys.argv[1])


