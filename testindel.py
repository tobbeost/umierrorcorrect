#!/usr/bin/env python3
from collections import Counter

def test(consensus_read,g):
    fsizes=[0,1,2,3,4,5,7,10,20,30]
    cons={}
    if consensus_read.indel_read==0:
        pos=consensus_read.start_pos
        count=consensus_read.count
        for j,base in enumerate(consensus_read.seq):
            if pos not in cons:
                cons[pos]={}
            for fsize in fsizes:
                if fsize==0:
                    if fsize not in cons[pos]:
                        cons[pos][fsize]=Counter()
                    cons[pos][fsize][base]+=count
                elif count>=fsize:
                    if fsize not in cons[pos]:
                        cons[pos][fsize]=Counter()
                    cons[pos][fsize][base]+=1
                    if fsize==10:
                        g.write('{} {} {} {}\n'.format(fsize,pos,j,base))
            pos+=1
    else:
        i=0
        pos=consensus_read.start_pos
        count=consensus_read.count
        cigar=consensus_read.cigarstring
        for base in consensus_read.seq:
            c=cigar[i]
            if pos not in cons:
                cons[pos]={}
            if c=='0':
                for fsize in fsizes:
                    if fsize==0:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize][base]+=count
                    elif count>=fsize:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize][base]+=1
                        if fsize==10:
                           g.write('{} {} {} {} {}\n'.format(fsize,pos,i,c,base))
                pos+=1
                i+=1
            elif c=='1':
                for fsize in fsizes:
                    if fsize==0:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize]['I']+=count
                    elif count>=fsize:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize]['I']+=1
                        if fsize==10:
                           g.write('{} {} {} {} {}\n'.format(fsize,pos,i,c,'I'))
                i+=1
            elif c=='2':
                for fsize in fsizes:
                    if fsize==0:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize]['D']+=count
                        
                    elif count>=fsize:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize]['D']+=1
                        if fsize==10:
                           g.write('{} {} {} {} {}\n'.format(fsize,pos,i,c,'D'))
                deletion=False        
                if cigar[i+1]=='2':
                    deletion=True
                    while deletion:
                        i+=1
                        c=cigar[i]
                        pos+=1
                        if pos not in cons:
                            cons[pos]={}
                        for fsize in fsizes:
                            if fsize==0:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize]=Counter()
                                cons[pos][fsize]['D']+=count
    
                            elif count>=fsize:
                                if fsize not in cons[pos]:
                                    cons[pos][fsize]=Counter()
                                cons[pos][fsize]['D']+=1
                            if fsize==10:
                               g.write('{} {} {} {} {}\n'.format(fsize,pos,i,c,'D'))    
                        if cigar[i+1]=='2':
                            deletion=True
                        else:
                            deletion=False
                pos+=1
                if pos not in cons:
                    cons[pos]={}
                for fsize in fsizes:
                    if fsize==0:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize][base]+=count

                    elif count>=fsize:
                        if fsize not in cons[pos]:
                            cons[pos][fsize]=Counter()
                        cons[pos][fsize][base]+=1
                        if fsize==10:
                           g.write('{} {} {} {} {}\n'.format(fsize,pos,i,c,base))
                pos+=1
                i+=1

def main():
    import pickle
    f=open('testnew5/cons.pickle','rb')
    cons_seq=pickle.load(f)
    for read in cons_seq.values():
        c=read.cigarstring
        if 'DD' in c:
            print(read.seq,read.name)
        if read.name=='Consensus_read_9_AGCGCCGAGCGT_Count=65':
            r2=read
        if read.name=='Consensus_read_9_ATAGTGGGTTGA_Count=52':
            r=read
        if read.name=='Consensus_read_9_TTCCTTAACCTC_Count=65':
            r0=read
        if read.name=='Consensus_read_9_AAGAAACAAAAT_Count=32':
            r3=read
    with open('indel1.out','w') as g:
        test(r2,g)
    with open('indel2.out','w') as g:
        test(r,g)
    with open('indel3.out','w') as g:
        test(r3,g)
    with open('indel4.out','w') as g:
        test(r0,g)

if __name__=='__main__':
    main()
