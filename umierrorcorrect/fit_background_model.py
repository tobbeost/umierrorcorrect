#!/usr/bin/env python3
import numpy as np
from numpy import inf
from numpy import nan
from scipy.optimize import fmin
from scipy.stats import beta
from scipy.special import beta as B
from scipy.special import comb
import argparse
import sys

def parseArgs():
    '''Function for parsing arguments'''
    parser = argparse.ArgumentParser(description="Pipeline for analyzing barcoded amplicon \
                                                  sequencing data with Unique molecular \
                                                  identifiers (UMI)")
    parser.add_argument('-cons', '--cons_file', dest='cons_file', help='Path to cons file, for fitting parameters of the bgmodel')
    parser.add_argument('-nonbgposfile', '--non-background-positions', dest='nonbgposfile',
                        help='Path to file with non-background positions')
    parser.add_argument('-out', '--out_file',dest='out_file',help="name of output file, default = %(default)s]",default="bgmodel.params") 
    parser.add_argument('-f','--fsize',dest='fsize', help='Family size cutoff (consensus cutoff) for variant calling. [default = %(default)s]', default=3)
    args = parser.parse_args(sys.argv[1:])
    return(args)

def parse_cons_file(filename,fsize=3):
    n1=[]
    f1=[]
    c1=[]
    posx=[]
    data=[]
    with open(filename) as f:
        for line in f:
            if not line.startswith('Sample Name'):
                line=line.rstrip('\n')
                parts=line.split('\t')
                pos=parts[1]+':'+parts[2]
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
                            posx.append(pos)
                            data.append(line)
                #print(name)
                #print(famsize)
    return(f1,n1,c1,posx,data)

def betaNLL(params,*args):
    a,b = params
    data = np.array(args[0])
    pdf=beta.pdf(data,a,b,loc=0,scale=1)
    lg=np.log(pdf)
    #lg=np.where(lg==-np.inf,0,lg)
    mask = np.isfinite(lg)
    nll = -lg[mask].sum()
    nll=-1*np.sum(lg)
    return(nll)


def get_beta_parameters(data):
    m=np.mean(data)
    v=np.var(data)
    a0=m*(m * (1-m) / v-1 )
    b0=(1-m)*(m * (1-m) / v-1 )
    result=fmin(betaNLL,[a0,b0],args=(data,))
    return(result)

def run_fit_bgmodel(args):
    spikepositions=[178952085,55599321,7577558,7577547,7577538,7577120]
    if args.nonbgposfile:
        nonbgpos=[]
        with open(args.nonbgposfile) as f:
            for line in f:
                line=line.rstrip()
                nonbgpos.append(line)
    else:
        nonbgpos=spikepositions
    if not args.cons_file:
        args.cons_file=glob.glob(args.output_path+'/*cons.tsv')[0]
    args.fsize=int(args.fsize)
    f1,n1,a1,pos,data=parse_cons_file(args.cons_file,args.fsize)
    f1 = np.array(f1)
    n1 = np.array(n1)
    a1 = np.array(a1)
    pos = np.array(pos)
    data = np.array(data)
    result=get_beta_parameters(f1[np.isin(pos,nonbgpos)!=True])
    #a=prob_bb(n1,a1,result[0],result[1])
    print(pos,nonbgpos,np.isin(pos,nonbgpos))
    with open(args.out_file,'w') as g:
        g.write('{}\n'.format(result[0]))
        g.write('{}\n'.format(result[1]))
    #a[a==inf]=1e-10
    #a[np.isnan(a)]=1e-10
    #Q = -10*np.log10(a)
    #data=np.array(data)
    #plot_histogram(Q,args.output_path+'/'+args.sample_name+'.histogram.png')
    #if args.vc_method.lower()=='bbmodel':
    #    rout=data[Q >= float(args.qvalue_threshold)]
    #    Qsig=Q[Q >= float(args.qvalue_threshold)]
    #else:
    #    rout=data[a1 >= float(args.count_cutoff)]
    #    Qsig=Q[a1 >= float(args.count_cutoff)]
    #outfilename=args.output_path+'/'+args.sample_name+'2.vcf'
    #write_vcf(outfilename,rout,Qsig,args.reference_file)

if __name__=='__main__':
    args=parseArgs()
    run_fit_bgmodel(args)
