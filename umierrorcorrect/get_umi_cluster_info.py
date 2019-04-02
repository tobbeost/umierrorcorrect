#!/usr/bin/env python3

import argparse
import sys
import logging

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)")
    parser.add_argument('-b', '--bam', dest='bam_file', help='Path to BAM-file')
    parser.add_argument('-d', '--edit_distance', dest='edit_distance_threshold',
                        help="Edit distance threshold for UMI clustering, [default = %(default)s]", default=1)
    parser.add_argument('-p', '--position_threshold', dest='position_threshold',
                        help='Position threshold for grouping by position [default = %(default)s]', default=10)
    parser.add_argument('-rid','--region_id', dest='region_id',help='region id')
    parser.add_argument('-bc','--barcode',dest='barcode', help='barcode sequence')
    parser.add_argument('-hist', dest='hist_file', help='hist file')
    args = parser.parse_args(sys.argv[1:])

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMI Error Correct')

    return(args)

def get_stat(stat_filename):
    with open(stat_filename) as f:
        regions=[]
        for line in f:
            line=line.rstrip()
            regionid, pos, name, cons, singles, *rest =line.split('\t')
            singles=int(singles.split(": ")[-1])
            regionid=int(regionid)
            regions.append((regionid,pos,singles,name))

    print(regions)

def main(args):
    regions=get_stat(args.hist_file)
    print(regions)
if __name__=='__main__':
    args=parseArgs()
    main(args)
