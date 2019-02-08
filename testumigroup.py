#!/usr/bin/env python3
from src.group import readBam

def main():
    bed_file='/medstore/Illumina_Tobias/newbatch/newbed3.bed'
    bam_file='/medstore/Illumina_Tobias/newbatch/A01_3B7/output.sorted.bam'
    #group_method='automatic'
    position_threshold=10
    regions, ends = readBam(bam_file, position_threshold)
    print(regions)

if __name__=='__main__':
    main()
