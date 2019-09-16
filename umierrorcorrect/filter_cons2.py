#!/usr/bin/env python3
import sys
import argparse

def parseArgs():
    parser=argparse.ArgumentParser(description="Filter cons file by removing positions without annotation and with a raw sequencing depth lower than threshold.")
    parser.add_argument('-i', '--infile',dest='infile', help='Path to the input file, required', required=True)
    parser.add_argument('-d','--raw_depth', dest='raw_depth_cutoff',help='Raw depth cutoff, [default = %(default)s]',default='150')
    parser.add_argument('-f','--family_sizes', dest='family_sizes',help='Family sizes to include, separated by comma. default= %(default)s',default='0,1,2,3,4,5,7,10,20,30')
    parser.add_argument('-only_annot','--only_annotation',dest='only_annotation',action='store_true', help='Include this flag if only annotated positions should be written to the output file')
    parser.add_argument('-write_raw', dest='writeraw',action='store_true',help='include this flag if raw reads should be included')
    args=parser.parse_args(sys.argv[1:])
    return(args)


def filter_cons(filename,raw_depth_cutoff=150,fsizes='0,1,2,3,4,5,7,10,20,30', writeraw=False, only_annotation=False):
    outfilename=filename[:-5]+'_filtered.cons'
    fs=fsizes.split(',')
    with open(filename) as f, open(outfilename,'w') as g:
        header=f.readline()
        g.write(header)
        passdepth=False
        for line in f:
            parts=line.split('\t')
            if only_annotation:
                if parts[3] not in '':
                    pos=parts[2]
                    fsize=parts[13]
                    if fsize=='0':
                        depth=int(parts[12])
                        if depth >= raw_depth_cutoff:
                            if writeraw:
                                g.write(line)
                            passdepth=True
                        else:
                            passdepth=False
                    elif passdepth==True and fsize in fs:
                        g.write(line)
            else:
                pos=parts[2]
                fsize=parts[13]
                if fsize=='0':
                    depth=int(parts[12])
                    if depth >= raw_depth_cutoff:
                        if writeraw:
                            g.write(line)
                        passdepth=True
                    else:
                        passdepth=False
                elif passdepth==True and fsize in fs:
                    g.write(line)

if __name__=='__main__':
    args=parseArgs()
    args.raw_depth_cutoff=int(args.raw_depth_cutoff)
    filter_cons(args.infile,args.raw_depth_cutoff,args.family_sizes,args.writeraw,args.only_annotation)


