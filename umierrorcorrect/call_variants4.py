#!/usr/bin/env python3
import numpy as np
from numpy import inf
from numpy import nan
from scipy.optimize import fmin
from scipy.stats import beta
from scipy.special import beta as B
from scipy.special import comb
import sys
import matplotlib.pyplot as plt

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

def prob_bb(n,k,a,b):
    p = comb(n,k) * B(k+a, n-k+b) / B(a,b)
    return(p)

def parse_cons_file(filename,fsize=3):
    n1=[]
    f1=[]
    data=[]
    with open(filename) as f:
        f.readline()
        for line in f:
            line=line.rstrip('\n')
            parts=line.split('\t')
            name=parts[2]
            #print(name)
            #if name not in "":
            famsize=parts[-3]
            if int(famsize)==fsize:
                frac=float(parts[-2])
                alt=parts[-1]
                if frac > 0 and alt not in 'N':
                    cov=int(parts[-4])
                    f1.append(float(frac))
                    n1.append(int(cov))
                    data.append(line)
                #print(name)
                #print(famsize)
    return(f1,n1,data)

def write_vcf(vcffile,rout,Qsig,reference):
    with open(vcffile,'w') as g:
        g.write('##fileformat=VCFv4.2\n##reference='+reference+ \
                '\n##source=umierrorcorrectV0.1\n##INFO=<ID=DP,Number=1,\
                Type=Integer,Description="Total Depth>\n')

def plot_histogram(hist,plot_filename):
    
    num_bins=100
    n,bins,patches=plt.hist(hist,num_bins,facecolor='dodgerblue',alpha=0.5)
    plt.xlabel('Q-score')
    plt.ylabel('Frequency')
    plt.title('Histogram of Q-scores')
    plt.box(False)
    plt.xlim(0,100)
    plt.savefig(plot_filename)


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
    plot_histogram(Q,'histogram.png')
    rout=data[Q>=cutoff]
    Qsig=Q[Q>=cutoff]
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
    main(sys.argv[1],int(sys.argv[2]),float(sys.argv[3]))
