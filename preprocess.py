#!/usr/bin/env python3
import sys
import gzip
from handle_sequences2 import read_fastq
def trim_barcode(sequence, barcode_length, spacer_length):
    barcode=sequence[:barcode_length]
    #spacer=sequence[barcode_length:barcode_length+spacer_length]
    rest=sequence[barcode_length+spacer_length:]
    return((barcode,rest))

def preprocess(infilename,outfilename):
    with gzip.open(infilename,'rb') as f, gzip.open(outfilename,'wb') as g:
        for name,seq,qualname,qual in read_fastq(f):
            #print (name,seq)
            barcode,rest=trim_barcode(seq,12,13)
            #g.write(name+':'+barcode+'\n'+rest+'\n'+qualname+'\n'+qual[12+11:]+'\n')
            newname=b':'.join([name,barcode])
            g.write(b'\n'.join([newname,rest,qualname,qual[12+11:]])+b'\n')

if __name__=='__main__':
    preprocess(sys.argv[1],sys.argv[2])
