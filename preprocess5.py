#!/usr/bin/env python3
import sys
import gzip
from handle_sequences2 import read_fastq, read_fastq_paired_end
import argparse
import os
import logging


def parseArgs():
    parser=argparse.ArgumentParser(description="Pipeline for analyzing  barcoded amplicon sequencing data with Unique molecular identifiers (UMI)3")
    parser.add_argument('-o', '--output_path',dest='output_path', help='Path to the output directory, required', required=True)
    parser.add_argument('-r1','--read1', dest='read1',help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2','--read2', dest='read2',help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-ul','--umi_length',dest='umi_length',help='Length of UMI sequence (number of nucleotides). The UMI is assumed to be located at the start of the read. Required',required=True)
    parser.add_argument('-sl','--spacer_length',dest='spacer_length',help='Length of spacer (The number of nucleotides between the UMI and the beginning of the read). The UMI + spacer will be trimmed off, and the spacer will be discarded. Default=%(default)s',default='0')
    
    parser.add_argument('-mode','--mode', dest='mode', help="Name of library prep, Either 'single' or 'paired', for single end or paired end data respectively, [default = %(default)s]",default="paired")
    parser.add_argument('-dual','--dual_index',dest='dual_index',help='Include this flag if dual indices are used (UMIs both on R1 and R2)',action='store_true')
    parser.add_argument('-reverse','--reverse_index',dest='reverse_index',help="Include this flag if a single index (UMI) is used, but the UMI is located on R2 (reverse read). Default is UMI on R1.", action='store_true')
    parser.add_argument('-tmpdir','--tmp_dir',dest='tmpdir',help="temp directory where the chunck files are written and then removed. Should be the scratch directory on the node. Only used if num_threads > 1. [default = %(default)s]", default="/tmp/")
    parser.add_argument('-cs','--chunk_size',dest='chunksize', help="Chunk size for reading the fastq files in chunks. Only used if num_threads > 1. [default = %(default)i]", default=25000)
    parser.add_argument('-t','--num_threads',dest='num_threads',help='Number of threads to run the program on. Default=%(default)s',default='1')
    args=parser.parse_args(sys.argv[1:])
    return(args)

def check_output_directory(outdir):
    if os.path.isdir(outdir):
        return(outdir)
    else:
        os.mkdir(outdir)
        return(outdir)

def generate_random_dir(tmpdir):
    import datetime
    newtmpdir=tmpdir+'/r'+datetime.datetime.now().strftime("%y%m%d_%H%M%S")+'/'
    newtmpdir=check_output_directory(newtmpdir)
    return(newtmpdir)

def get_sample_name(read1,mode):
    if mode=='single':
        samplename=read1.split('/')[-1].rstrip('fastq').rstrip('fastq.gz')
    elif mode=='paired':
        samplename=read1.split('/')[-1].rstrip('fastq').rstrip('fastq.gz').rstrip('_R012')
    return(samplename)

def chunks_paired(read1,read2,chunksize,tmpdir):
    fid = 1
    name_list = []
    with gzip.open(read1,'rb') as infile1:#, gzip.open(read2,'rb') as infile2:
        chunkname = tmpdir + '/' +  'chunk%d' % fid
        f1 = gzip.open(chunkname+'_1.fastq.gz', 'wb')
        #f2 = gzip.open(chunkname+'_2.fastq.gz', 'wb')
        for i, a in enumerate(infile1):
            f1.write(a)
            #f2.write(b)
            
            if not (i+1) % (chunksize*4) and not i==0:
                f1.close()
                #f2.close()
                name_list.append(chunkname+'_1.fastq.gz')
                fid += 1
                chunkname = tmpdir + '/chunk%d' % fid
                f1 = gzip.open(chunkname+'_1.fastq.gz', 'wb')
                #f2 = gzip.open(chunkname+'_2.fastq.gz', 'wb')
        #name_list.append((chunkname+'_1.fastq.gz',chunkname+'_2.fastq.gz'))
        name_list.append(chunkname+'_1.fastq.gz')
        f1.close()
    return name_list



def trim_barcode(sequence, barcode_length, spacer_length):
    barcode=sequence[:barcode_length]
    #spacer=sequence[barcode_length:barcode_length+spacer_length]
    rest=sequence[barcode_length+spacer_length:]
    return((barcode,rest))

def preprocess_se(infilename,outfilename,barcode_length,spacer_length):
    try:
        barcode_length=int(barcode_length)
    except ValueError as e:
        logger.warning(e+" Barcode length needs to be an integer")
        sys.exit(1)
    try:
        spacer_length=int(spacer_length)
    except ValueError as e:
        logger.warning(e+" Spacer length needs to be an integer")
        sys.exit(1) 
    
    with gzip.open(infilename,'rb') as f, gzip.open(outfilename,'wb') as g:
        read_start=barcode_length+spacer_length
        for name,seq,qual in read_fastq(f):
            barcode=seq[:barcode_length]
           #g.write(name+':'+barcode+'\n'+rest+'\n'+qualname+'\n'+qual[12+11:]+'\n')
            parts=name.split()
            newname=b':'.join([parts[0],barcode])+b' '+parts[-1]
            g.write(b'\n'.join([newname,seq[read_start:],b'+',qual[read_start:]])+b'\n')

def preprocess_pe(r1file,r2file,outfile1,outfile2,barcode_length,spacer_length,dual_index):
    try:
        barcode_length=int(barcode_length)
    except ValueError as e:
        logger.warning(e+" Barcode length needs to be an integer")
        sys.exit(1)
    try:
        spacer_length=int(spacer_length)
    except ValueError as e:
        logger.warning(e+" Spacer length needs to be an integer")
        sys.exit(1)
    read_start=barcode_length+spacer_length
    with gzip.open(r1file,'rb') as f1, gzip.open(r2file,'rb') as f2, gzip.open(outfile1,'wb') as g1, gzip.open(outfile2,'wb') as g2:
        for name1,seq1,qual1,name2,seq2,qual2 in read_fastq_paired_end(f1,f2):
            if dual_index==True:
                barcode=seq1[:barcode_length]+seq2[:barcode_length]
            else:
                barcode=seq1[:barcode_length]
            parts1=name1.split()
            parts2=name2.split()
            newname1=b':'.join([parts1[0],barcode])+b' '+parts1[-1]
            newname2=b':'.join([name2,barcode])+b' '+parts2[-1]
            g1.write(b'\n'.join([newname1,seq1[read_start:],b'+',qual1[read_start:]])+b'\n')
            if dual_index==True:
                g2.write(b'\n'.join([newname2,seq2[read_start:],b'+',qual2[read_start:]])+b'\n')
            else:
                g2.write(b'\n'.join([newname2,seq2,b'+',qual2])+b'\n')


    #if dual_index==False:
    #    with gzip.open(r1file,'rb') as f, gzip.open(outfile1,'wb') as g:
    #        for name,seq,qual in read_fastq(f):
    #            barcode=seq[:barcode_length]
    #            read_id=name.split()[0]
    #            if read_id not in fastqdict:
    #                fastqdict[read_id]=barcode
    #            #else:
    #            #   print('duplicate')
    #            newname=b':'.join([name,barcode])
    #            g.write(b'\n'.join([newname,seq[read_start:],b'+',qual[read_start:]])+b'\n') 
    #    with gzip.open(r2file,'rb') as f, gzip.open(outfile2,'wb') as g:
    #        for name,seq,qual in read_fastq(f):
    #            read_id=name.split()[0]
    #            if read_id in fastqdict:
    #                barcode=fastqdict[read_id]
    #                newname=b':'.join([name,barcode])
    #                g.write(b'\n'.join([newname,seq,b'+',qual])+b'\n')
    #            #else:
    #            #    print("missing")
    #
    #elif dual_index==True:
    #    seqlist=[]
    #    fastqdict2={}
    #    with gzip.open(r1file,'rb') as f:
    #        for name,seq,qual in read_fastq(f):
    #            barcode=seq[:barcode_length]
    #            read_id=name.split()[0]
    #            if read_id not in fastqdict:
    #                fastqdict[read_id]=barcode
    #            #else:
    #            #    print('duplicate')
    #            seqlist.append((name,seq[read_start:],qual[read_start:]))
    #    with gzip.open(r2file,'rb') as f, gzip.open(outfile2,'wb') as g:
    #        for name,seq,qual in read_fastq(f):
    #            read_id=name.split()[0]
    #            if read_id in fastqdict:
    #                barcode=seq[:barcode_length]
    #                newbarcode=fastqdict[read_id]+barcode
    #                if read_id not in fastqdict2:
    #                    fastqdict2[read_id]=barcode
    #                newname=b':'.join([name,newbarcode])
    #                g.write(b'\n'.join([newname,seq[read_start:],b'+',qual[read_start:]])+b'\n')
    #            #else:
    #            #    print(missing)
    #    with gzip.open(outfile1,'wb') as g:
    #        for name,seq,qual in seqlist:
    #            read_id=name.split()[0]
    #            newbarcode=fastqdict[read_id]+fastqdict2[read_id]
    #            newname=b':'.join([name,newbarcode])
    #            g.write(b'\n'.join([newname,seq,b'+',qual])+b'\n')



def main(args):
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S',level=logging.DEBUG)
    logging.info('Starting UMI Error Correct')
    args.output_path=check_output_directory(args.output_path)
    newtmpdir=generate_random_dir(args.tmpdir)
    args.chunksize=int(args.chunksize)
    #if int(args.num_threads)>1:
        #if args.mode=='paired':
            #filelist=chunks_paired(args.read1,args.read2,args.chunksize,newtmpdir)
            #print(filelist)
    logging.info('Writing output files to {}'.format(args.output_path))
    if args.mode=='paired':
        if not args.read2:
            logging.warning("R1 and R2 Fastq files are required for mode 'paired', exiting.")
            sys.exit(1)
    samplename=get_sample_name(args.read1,args.mode)
    if args.mode=='single':
        outfilename=args.output_path+'/'+samplename+'.fastq.gz'
        preprocess_se(args.read1,outfilename,args.umi_length,args.spacer_length)
    else:
        if args.reverse_index==True:
            #switch forward and reverse read
            r1file=args.read2
            r2file=args.read1
            outfile1=args.output_path+'/'+samplename+'_R2.fastq.gz'
            outfile2=args.output_path+'/'+samplename+'_R1.fastq.gz'
        else:
            r1file=args.read1
            r2file=args.read2
            outfile1=args.output_path+'/'+samplename+'_R1.fastq.gz'
            outfile2=args.output_path+'/'+samplename+'_R2.fastq.gz'
        preprocess_pe(r1file,r2file,outfile1,outfile2,args.umi_length,args.spacer_length,args.dual_index)
if __name__=='__main__':
    args=parseArgs()
    main(args)
