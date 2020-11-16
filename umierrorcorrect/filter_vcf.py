#!/usr/bin/env python3
import sys
import argparse

def parseArgs():
    parser=argparse.ArgumentParser(description="Filter vcf file by removing positions llower than threshold.")
    parser.add_argument('-i', '--infile',dest='infile', help='Path to the input file, required', required=True)
    parser.add_argument('-Q', '--qscore_cutoff', dest='qvalue_threshold',
                                help='Qscore threshold (Minimum variant significance score) for Variant calling [default = %(default)s]',
                                                        default=20)
    parser.add_argument('-d','--depth', dest='depth_cutoff',help='Depth cutoff, [default = %(default)s]',default='5')
    parser.add_argument('-filter', dest='filter', help='filter that should be included')
    args=parser.parse_args(sys.argv[1:])
    return(args)


def filter_cons(filename,filter,depth_cutoff=5,qvalue_threshold=20):
    outfilename=filename[:-4]+'_filtered.vcf'

    with open(filename) as f, open(outfilename,'w') as g:
        for line in f:
            if line.startswith('#'):
                g.write(line)
            else:
                line=line.rstrip()
                parts=line.split('\t')
                Q=float(parts[5])
                f=parts[6]
                d=int(parts[-1])
                if filter==None:
                    if Q >= qvalue_threshold and d >= depth_cutoff:
                        g.write(line+'\n')
                else:
                    if f==filter and Q >= qvalue_threshold and d >= depth_cutoff:
                        g.write(line+'\n')


if __name__=='__main__':
    args=parseArgs()
    args.depth_cutoff=int(args.depth_cutoff)
    args.qvalue_threshold=int(args.qvalue_threshold)
    filter_cons(args.infile,args.filter,args.depth_cutoff,args.qvalue_threshold)
