#!/usr/bin/env python3
#import numpy as np
#from numpy import inf
#from numpy import nan
#from scipy.optimize import fmin
#from scipy.stats import beta
#from scipy.special import beta as B
#from scipy.special import comb
import argparse
import glob
import sys
import logging
import subprocess
#import matplotlib.pyplot as plt

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)")
    parser.add_argument('-o',  '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-cons', '--cons_file', dest='cons_file', help='Path to cons file')
    parser.add_argument('-s', '--sample_name', dest='sample_name',
                        help='Sample name that will be used as base name for the output files. \
                        If excluded the sample name will be extracted from the BAM file.')
    parser.add_argument('-r', '--reference', dest='reference_file',
                            help='Path to the reference sequence in Fasta format, Used for annotation', default='reference')
    parser.add_argument('-f','--fsize',dest='fsize', help='Family size cutoff (consensus cutoff) for variant calling. [default = %(default)s]', default=3)
    parser.add_argument('-method','--vc-method',dest='vc_method',
                        help="Variant calling method, Either 'count' or 'bbmodel'. [default = %(default)s]", default='count')
    parser.add_argument('-count', '--count_cutoff', dest='count_cutoff',
                        help="Consensus read count cutoff (minimum variant allele depth) for calling a variant if method=count [default = %(default)s]", default=5)
    parser.add_argument('-Q', '--qscore_cutoff', dest='qvalue_threshold',
                        help='Qscore threshold (Minimum variant significance score) for Variant calling, only if method=bbmodel [default = %(default)s]', default=20)
    args = parser.parse_args(sys.argv[1:])

    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMI Error Correct')

    return(args)


def runR(consfilename, outfilename, fsize):
    command=['Rscript','/home/xsteto/umierrorcorrect/umierrorcorrect/call_variants2.R',consfilename,outfilename,str(fsize)]
    p1 = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.communicate()

def parse_R_output(filename):
    Rout=[]
    with open(filename) as f:
        f.readline()
        for line in f:
            line=line.rstrip('\n')
            parts=line.split('\t')
            Rout.append(line)
    return(Rout)

def write_vcf(vcffile,rout,reference):
    with open(vcffile,'w') as g:
        g.write('##fileformat=VCFv4.2\n##reference='+reference+ \
                '\n##source=umierrorcorrectV0.1\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'+\
                '##INFO=<ID=AD,Number=1,Type=Float,Description="Alternative Allele Depth">\n'+\
                '##INFO=<ID=AF,Number=A,Type=Integer,Description="Alternative Allele Frequency">\n'+\
                '##FILTER=<ID=a5,Description="Alternative Allele Depth below 5">\n'+\
                '##FILTER=<ID=q10,Description="Variant quality below 10">\n'+\
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Allele Depth">\n')
        g.write('\t'.join(['#CHROM',  'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])+'\n')

        for r in rout:
            parts=r.split('\t')
            q=parts[-1]
            if q=='NA':
                q='.'
            vcffilter=[]
            if not q=='.' and float(q) > 0.05:
                vcffilter.append('q10')
            if int(parts[-5]) < 5:
                vcffilter.append('a5')
            if len(vcffilter)==0:
                vcffilter='PASS'
            else:
                vcffilter=','.join(vcffilter)
            newline='\t'.join([parts[1], parts[2], '.', parts[4], parts[-3], str(q), vcffilter, 'DP={};AD={};AF={}'.format(parts[-7],parts[-5],parts[-4]),'DP',parts[-5]])
            g.write(newline+'\n')


def get_sample_name(cons_name):
    '''Get the sample name as the basename of the input files.'''
    sample_name=cons_name.split('/')[-1]
    sample_name = sample_name.replace('.cons','')
    return(sample_name)

def run_call_variants(args):
    if not args.cons_file:
        args.cons_file=glob.glob(args.output_path+'/*.cons')[0]
    if not args.sample_name:
        args.sample_name=get_sample_name(args.cons_file)
    runR(args.cons_file,args.output_path+'/output.txt',args.fsize)
    outfilename=args.output_path+'/'+args.sample_name+'.vcf'
    Rout= parse_R_output(args.output_path+'/output.txt')
    write_vcf(outfilename,Rout,args.reference_file)

def main(filename,fsize,cutoff):
    f1,n1,data=parse_cons_file(filename,fsize)
    print(f1)
    f1 = np.array(f1)
    n1 = np.array(n1)
    a1 = f1*n1
    result=get_beta_parameters(f1)
    a=prob_bb(n1,a1,result[0],result[1])
    a[a==inf]=1e-10
    a[np.isnan(a)]=1e-10
    Q = -10*np.log10(a)
    data=np.array(data)
    #plot_histogram(Q,'histogram.png')
    rout=data[Q>=float(args.qvalue_threshold)]
    Qsig=Q[Q>=float(args.qvalue_threshold)]
    outfilename=filename.rstrip('.cons')+'.vcf'
    write_vcf(outfilename,rout,Qsig,'reference')
    for r,q in zip(rout,Qsig):
        print(r+'\t'+str(q))

    #with open('fraction.txt') as f:
    #    data=[]
    #    for line in f:
    #        line=line.rstrip()
    #        data.append(float(line))
    #result=beta.fit(f1,floc=0)
	#
    #result=get_beta_parameters(f1)
    #print(result)

if __name__=='__main__':
    args=parseArgs()
    run_call_variants(args)
    #main(sys.argv[1],int(sys.argv[2]),float(sys.argv[3]))
