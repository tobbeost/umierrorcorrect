#!/usr/bin/env python3
import sys

def main(filename1,filename2):
    #filename1='umi_no_substr3.txt'
    #filename2='umi_substr3.txt'
    with open(filename1) as f1, open(filename2) as f2:
        lines1=f1.readlines()
        lines2=f2.readlines()
        print (len(lines1),len(lines2))
        for l1,l2 in zip(lines1,lines2):
            if not l1==l2:
                print(l1,l2)
if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])
