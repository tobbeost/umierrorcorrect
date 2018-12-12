#!/bin/bash -l
#$ -cwd
#$ -o run_umierrorcorrect_181205.stdout
#$ -e run_umierrorcorrect_181205.stderr
#$ -q production.q@hoth.medair.lcl
#$ -pe mpi 40
#S -l excl=1
module load samtools
source activate debarcer
#mkdir /tmp/tobiastmp/
#cp GMS_solid/C6280_S7_R* /tmp/tobiastmp/
#time ./preprocess5.py -r1 /tmp/tobiastmp/C6280_S7_R1_001.fastq.gz -r2 /tmp/tobiastmp/C6280_S7_R2_001.fastq.gz -o /tmp/tobiastmp/outtest5 -ul 3 -sl 2 -dual
time bwa mem -t 40 /medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta /tmp/tobiastmp/outtest5/C6280_S7_R1.fastq.gz /tmp/tobiastmp/outtest5/C6280_S7_R2.fastq.gz | samtools view -Sb -@40 -o /tmp/tobiastmp/output.bam
samtools sort -@40 /tmp/tobiastmp/output.bam -o /tmp/tobiastmp/output.sorted.bam
samtools index -@40 /tmp/tobiastmp/output.sorted.bam
cp -r /tmp/tobiastmp/output.bam* .

