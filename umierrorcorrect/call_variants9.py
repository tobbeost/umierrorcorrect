#!/usr/bin/env python3
import numpy as np
from numpy import inf
from numpy import nan
from scipy.optimize import fmin
from scipy.stats import beta
from scipy.special import beta as B
from scipy.special import comb
import argparse
import glob
import sys
import logging
import matplotlib.pyplot as plt

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)")
    parser.add_argument('-o',  '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-cons', '--cons_file', dest='cons_file', help='Path to cons file')
    parser.add_argument('-r', '--reference', dest='reference_file',
                        help='Path to the reference sequence in Fasta format, Used for annotation, required', default='reference')
    parser.add_argument('-s', '--sample_name', dest='sample_name',
                        help='Sample name that will be used as base name for the output files. \
                        If excluded the sample name will be extracted from the BAM file.')
    parser.add_argument('-p','--params_file',dest='params_file',help='Params file')
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


def prob_bb(n,k,a,b):
    p = comb(n,k) * B(k+a, n-k+b) / B(a,b)
    return(p)

def parse_cons_file(filename,fsize=3):
    n1=[]
    f1=[]
    c1=[]
    data=[]
    with open(filename) as f:
        f.readline()
        for line in f:
            line=line.rstrip('\n')
            parts=line.split('\t')
            pos=parts[2]
            name=parts[3]
            #print(name)
            if name not in "":
                famsize=parts[-4]
                if int(famsize)==fsize:
                    frac=float(parts[-2])
                    alt=parts[-1]
                    count=parts[-3]
                    if frac > 0 and alt not in 'N':
                        cov=int(parts[-5])
                        f1.append(float(frac))
                        n1.append(int(cov))
                        c1.append(int(count))
                        data.append(line)
                #print(name)
                #print(famsize)
    return(f1,n1,c1,data)

def write_vcf(vcffile,rout,Qsig,reference):
    with open(vcffile,'w') as g:
        g.write('##fileformat=VCFv4.2\n##reference='+reference+ \
                '\n##source=umierrorcorrectV0.1\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'+\
                '##INFO=<ID=AD,Number=1,Type=Float,Description="Alternative Allele Depth">\n'+\
                '##INFO=<ID=AF,Number=A,Type=Integer,Description="Alternative Allele Frequency">\n'+\
                '##FILTER=<ID=a5,Description="Alternative Allele Depth below 5">\n'+\
                '##FILTER=<ID=q10,Description="Variant quality below 10">\n'+\
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Allele Depth">\n')
        g.write('\t'.join(['#CHROM',  'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])+'\n')

        for r,q in zip(rout,Qsig):
            if q=='NA':
                q='.'
            parts=r.split('\t')
            vcffilter=[]
            if not q=='.' and float(q) < 10:
                vcffilter.append('q10')
            if int(parts[-3]) < 5:
                vcffilter.append('a5')
            if len(vcffilter)==0:
                vcffilter='PASS'
            else:
                vcffilter=','.join(vcffilter)
            newline='\t'.join([parts[1], parts[2], '.', parts[4], parts[-1], str(q), vcffilter, 'DP={};AD={};AF={}'.format(parts[-5],parts[-3],parts[-2]),'DP',parts[-3]])
            g.write(newline+'\n')

def plot_histogram(hist,plot_filename):
    
    num_bins=100
    n,bins,patches=plt.hist(hist,num_bins,facecolor='dodgerblue',alpha=0.5)
    plt.xlabel('Q-score')
    plt.ylabel('Frequency')
    plt.title('Histogram of Q-scores')
    plt.box(False)
    plt.xlim(0,140)
    plt.savefig(plot_filename)

def get_sample_name(cons_name):
    '''Get the sample name as the basename of the input files.'''
    sample_name=cons_name.split('/')[-1]
    sample_name = sample_name.replace('_cons.tsv','')
    return(sample_name)

def run_call_variants(args):
    spikepositions=[178952085,55599321,7577558,7577547,7577538,7577120]
    if not args.cons_file:
        args.cons_file=glob.glob(args.output_path+'/*cons.tsv')[0]
    if not args.sample_name:
        args.sample_name=get_sample_name(args.cons_file)
    args.fsize = int(args.fsize)
    f1,n1,a1,data=parse_cons_file(args.cons_file,args.fsize)
    f1 = np.array(f1)
    n1 = np.array(n1)
    a1 = np.array(a1)
    data = np.array(data)
    if args.vc_method.lower()=='count':
        rout=data[a1 >= float(args.count_cutoff)]
        Qsig=['NA']*len(rout)
    #result=get_beta_parameters(f1[np.isin(pos,spikepositions)!=True])
    params=[]
    if args.params_file:
        with open(args.params_file) as f:
            for line in f:
                line=line.rstrip()
                params.append(float(line))
    else:
        params=[2.168215069116764,3531.588541594945]
    a=prob_bb(n1,a1,params[0],params[1])
    print(params[0])
    print(params[1])
    a[a==inf]=1e-10
    a[np.isnan(a)]=1e-10
    a[a==0]=1e-10
    Q = -10*np.log10(a)
    data=np.array(data)
    #plot_histogram(Q,args.output_path+'/'+args.sample_name+'.histogram.png')
    if args.vc_method.lower()=='bbmodel':
        rout=data[Q >= float(args.qvalue_threshold)]
        Qsig=Q[Q >= float(args.qvalue_threshold)]
    else:
        rout=data[a1 >= float(args.count_cutoff)]
        Qsig=Q[a1 >= float(args.count_cutoff)]
    outfilename=args.output_path+'/'+args.sample_name+'.vcf'
    write_vcf(outfilename,rout,Qsig,args.reference_file)

def main(filename,fsize,cutoff):
    spikepositions=[178952085,55599321,7577558,7577547,7577538,7577120]
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
    rout=data[Q>=float(args.qvalue_threshold)]
    Qsig=Q[Q>=float(args.qvalue_threshold)]
    outfilename=filename.rstrip('_cons.tsv')+'.vcf'
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
