#!/usr/bin/env python3
def read_fastq(infile):
    #reading one fastq-record at a time using a generator
    for line in infile:
        name=line.rstrip()
        seq=infile.readline().rstrip()
        infile.readline()
        qual=infile.readline().rstrip()
        yield(name,seq,qual)


def read_fastq_paired_end(r1file,r2file):
    for name1,seq1,qual1 in read_fastq(r1file): 
        name2,seq2,qual2=next(read_fastq(r2file))
        #if name1.split()[0] != name2.split()[0]:
        #    raise ValueError("Read names do not match for pair {},{}".format(name1,name2))
        yield name1,seq1,qual1,name2,seq2,qual2
    
