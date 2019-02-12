#!/usr/bin/env python3


def read_fastq(infile):
    # reading one fastq-record at a time using a generator
    name, seq, qual = None, None, None
    for line in infile:
        name = line.rstrip()
        seq = infile.readline().rstrip()
        infile.readline()
        # qualname=infile.readline().rstrip()
        qual = infile.readline().rstrip()
        yield(name, seq, qual)


def read_fastq_paired_end(r1file, r2file):
    for line1, line2 in zip(r1file, r2file):
        name1 = line1.rstrip()
        name2 = line2.rstrip()
        # assert name1.split()[0] == name2.split()[0]
        seq1 = r1file.readline().rstrip()
        seq2 = r2file.readline().rstrip()
        r1file.readline()
        r2file.readline()
        qual1 = r1file.readline().rstrip()
        qual2 = r2file.readline().rstrip()
        yield(name1, seq1, qual1, name2, seq2, qual2)
