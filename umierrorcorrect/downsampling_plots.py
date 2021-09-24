#!/usr/bin/env python3
import sys
import argparse
import logging
import glob
from umierrorcorrect.get_consensus_statistics import get_stat, downsample_reads_per_region, save_downsampled_table, region_cons_stat, plot_downsampling
from umierrorcorrect.version import __version__

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="UmiErrorCorrect v. {}. \
                                                  Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)".format(__version__))
    parser.add_argument('-o',  '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-c', '--cons_bam', dest='cons_bam_file', help='Path to the consensus BAM-file')
    parser.add_argument('-hist', '--hist_file', dest='hist_file',
                        help='Path to the .hist file')
    parser.add_argument('-s','--sample_name',dest='samplename', help='Sample name, if not provided it is extracted')
    parser.add_argument('-f','--fsize',dest='fsize', help='Family size cutoff (consensus cutoff) for downsampling. [default = %(default)s]', default='3')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMI Error Correct')
    return(args)

def run_downsampling(output_path, consensus_filename, stat_filename, fsize,samplename=None):
    logging.info('Getting consensus statistics')
    if not consensus_filename:
        consensus_filename=glob.glob(output_path + '/*_consensus_reads.bam')[0]
    if not samplename:
        samplename = consensus_filename.split('/')[-1].replace('_consensus_reads.bam','')
    if not stat_filename:
        stat_filename=output_path+'/'+samplename+'.hist'
    hist=get_stat(consensus_filename,stat_filename)
    fsizes=[1,2,3,4,5,7,10,20,30]
    tot_results=region_cons_stat('All','all_regions','',0,fsizes)
    for h in hist:
        tot_results.hist=tot_results.hist + h.hist
        tot_results.singletons+=h.singletons
        
    downsample_rates=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    tot=downsample_reads_per_region([tot_results],downsample_rates,fsizes,False)
    all_results=downsample_reads_per_region(hist, downsample_rates, fsizes, True)
    filename=output_path+'/'+samplename+'_downsampled_coverage.txt'
    save_downsampled_table(all_results,tot,filename)
    filename=output_path+'/'+samplename+'_downsampled_plot.png'
    plot_downsampling(tot,fsize,filename)

if __name__=='__main__':
    args=parseArgs()
    run_downsampling(args.output_path, args.cons_bam_file, args.hist_file, args.fsize, args.samplename)
