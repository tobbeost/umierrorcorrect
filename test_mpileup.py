#/usr/bin/env python3
import pysam
from collections import Counter
from get_cons_dict2 import get_reference_sequence
from get_regions_from_bed2 import get_annotation,read_bed,sort_regions,merge_regions

def run_mpileup(bamfilename,fsizes=[0,1,2,3,4,5,7,10,20,30]):
    cons={}
    with pysam.AlignmentFile(bamfilename,'rb') as f:
        alignment=f.pileup(max_depth=1000000)
        for pileupcolumn in alignment:
            pos=pileupcolumn.pos
            chrx=pileupcolumn.reference_name
            if chrx not in cons:
                cons[chrx]={}
            for read in pileupcolumn.pileups:
                count=int(read.alignment.qname.split('=')[-1])
                regionid=read.alignment.qname.split('_')[2]
                if regionid not in cons[chrx]:
                    cons[chrx][regionid]={}
                if pos not in cons[chrx][regionid]:
                    cons[chrx][regionid][pos]={}
                if not read.is_del and not read.indel:
                    base=read.alignment.query_sequence[read.query_position]
                elif read.indel < 0:
                    base='D'
                elif read.indel > 0:
                    base='I'
                for fsize in fsizes:
                    if fsize==0:
                        if fsize not in cons[chrx][regionid][pos]:
                            cons[chrx][regionid][pos][fsize]=Counter()
                        cons[chrx][regionid][pos][fsize][base]+=count
                    elif count>=fsize:
                        if fsize not in cons[chrx][regionid][pos]:
                            cons[chrx][regionid][pos][fsize]=Counter()
                        cons[chrx][regionid][pos][fsize][base]+=1 
    return(cons)


def calc_major_nonref_allele_frequency(cons,ref):
    tot=sum(cons.values())
    comp={key:cons[key] for key in cons if key != ref}
    allele=max(comp,key=comp.get)
    frac=1.0*(cons[allele]/tot)
    if frac > 0:
        return((allele,frac,tot))
    else:
        return(("",0,tot))

def write_consensus(f,cons,ref_seq,start,contig,annotation,only_target_regions):
    bases=['A','C','G','T','I','D','N']
    #print(list(cons.keys())[0],list(cons.keys())[-1],start,len(ref_seq))

    for pos in sorted(cons):
        annotation_pos=get_annotation(annotation,pos)
        if not (annotation_pos=="" and only_target_regions):
            #if len(ref_seq)<(pos-start+1):
            #    print("error",contig,start,ref_seq)
            pos=int(pos)
            start=int(start)
            refbase=ref_seq[pos-start]

            for fsize in cons[pos]:
                line=[]
                line.append(contig)
                line.append(str(pos))
                line.append(annotation_pos)
                line.append(refbase)
                if len(cons[pos][fsize])>1:
                    mna,freq,tot=calc_major_nonref_allele_frequency(cons[pos][fsize],refbase)
                else:
                    mna=''
                    freq=0
                    tot=sum(cons[pos][fsize].values())
                for base in bases:
                    if base in cons[pos][fsize]:
                        line.append(str(cons[pos][fsize][base]))
                    else:
                        line.append(str(0))
                line.append(str(tot))
                line.append(str(fsize))
                line.append(str(freq))
                line.append(mna)
                f.write('\t'.join(line)+'\n')


def main(bamfilename,bedfilename):
    cons=run_mpileup(bamfilename)
    regions=read_bed(bedfilename)
    regions=sort_regions(regions)
    regions=merge_regions(regions,0)
    with open(bamfilename[:-4]+'.cons','w') as g, pysam.FastaFile('/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta') as fasta:
        for contig in cons:
            for regionid in cons[contig]:
                start=min(list(cons[contig][regionid].keys()))
                end=max(list(cons[contig][regionid].keys()))+1
                print(contig,start,end)
                reference_sequence=get_reference_sequence(fasta,contig,start,end)
                #print(reference_sequence)
                if contig in regions:
                    annotation=regions[contig]
                else:
                    annotation=[]
                write_consensus(g,cons[contig][regionid],reference_sequence,start,contig,annotation,False)


if __name__=='__main__':
    import sys
    main(sys.argv[1],sys.argv[2])
