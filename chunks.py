#!/usr/bin/env python3
import sys
def chunks_paired(read1,read2,chunksize,tmpdir):
    fid = 1
    name_list = []
    with open(read1,'r') as infile1, open(read2,'r') as infile2:
        chunkname = tmpdir + '/' +  'chunk%d' % fid
        print(chunkname)
        print(chunksize)
        f1 = open(chunkname+'_1.fastq', 'w')
        f2 = open(chunkname+'_2.fastq', 'w')

        for i, (a ,b) in enumerate(zip(infile1,infile2)):
            f1.write(a)
            f2.write(b)

            if not (i+1) % (chunksize*4) and not i==0:
                print(i)
                f1.close()
                f2.close()
                name_list.append(chunkname+'_1.fastq')
                fid += 1
                chunkname = tmpdir + '/chunk%d' % fid
                f1 = open(chunkname+'_1.fastq', 'w')
                f2 = open(chunkname+'_2.fastq', 'w')
        name_list.append((chunkname+'_1.fastq',chunkname+'_2.fastq'))
        #name_list.append(chunkname+'_1.fastq.gz')
        f1.close()
        f2.close()
    return name_list


def main(r1,r2):
    filelist=chunks_paired(r1,r2,100000,'/tmp/tobiastmp/r123')
    print(filelist)

if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])
