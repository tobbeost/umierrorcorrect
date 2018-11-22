#!/usr/bin/env python3
def read_fastq(infile):
    #reading one fastq-record at a time using a generator
    name, seq, qualname, qual = None, None, None, None
    for line in infile:
        name=line.rstrip()
        seq=infile.readline().rstrip()
        infile.readline()
        #qualname=infile.readline().rstrip()
        qual=infile.readline().rstrip()
        yield(name,seq,qual)
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
def read_fastq_paired_end(r1file,r2file):
    while True:
        for line1,line2 in zip(r1file,r2file):
            name1=line1.rstrip()
            name2=line2.rstrip()
            break
        assert name1.split()[0] == name2.split()[0]
        for line1,line2 in zip(r1file,r2file):
            seq1=line1.rstrip()
            seq2=line2.rstrip()
            break
        for line1,line2 in zip(r1file,r2file):
            break
        for line1,line2 in zip(r1file,r2file):
            qual1=line1.rstrip()
            qual2=line2.rstrip()
            break    
        yield(name1,seq1,qual1,name2,seq2,qual2)
    
