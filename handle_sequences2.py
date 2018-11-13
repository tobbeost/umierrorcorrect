#!/usr/bin/env python3
def read_fastq(infile):
    #reading one fastq-record at a time using a generator
    name, seq, qualname, qual = None, None, None, None
    for line in infile:
        name=line.rstrip()
        seq=infile.readline().rstrip()
        qualname=infile.readline().rstrip()
        qual=infile.readline().rstrip()
        yield(name,seq,qualname,qual)

