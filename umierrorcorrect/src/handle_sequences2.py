#!/usr/bin/env python3


def read_fastq(infile):
    # reading one fastq-record at a time using a generator
    name, seq, qual = None, None, None
    while True:
        line = infile.readline()
        if not line:
            break
        name = line.rstrip()
        line = infile.readline()
        if not line:
            break
        seq = line.rstrip()
        line = infile.readline()
        if not line:
            break
        # qualname=infile.readline().rstrip()
        line = infile.readline()
        if not line:
            break
        qual = line.rstrip()
        yield(name, seq, qual)


def read_fastq_paired_end(r1file, r2file):
    while True:
        line1 = r1file.readline() 
        line2 = r2file.readline()
        if not line1 or not line2:
            break
        name1 = line1.rstrip()
        name2 = line2.rstrip()
        # assert name1.split()[0] == name2.split()[0]
        line1 = r1file.readline()
        line2 = r2file.readline()
        if not line1 or not line2:
            break
        seq1 = line1.rstrip()
        seq2 = line2.rstrip()
        line1 = r1file.readline()
        line2 = r2file.readline()
        if not line1 or not line2:
            break
        line1 = r1file.readline()
        line2 = r2file.readline()
        if not line1 or not line2:
            break
        qual1 = line1.rstrip()
        qual2 = line2.rstrip()
        yield(name1, seq1, qual1, name2, seq2, qual2)
