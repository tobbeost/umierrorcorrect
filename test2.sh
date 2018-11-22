#!/bin/bash -l
#$ -cwd
#$ -o run_umierrorcorrect_181122.stdout
#$ -e run_umierrorcorrect_181122.stderr
#$ -q production.q@hoth.medair.lcl
#$ -pe mpi 40
#S -l excl=1
source activate debarcer
mkdir /tmp/tobiastmp/
cp GMS_solid/C6280_S7_R* /tmp/tobiastmp/
time ./preprocess5.py -r1 /tmp/tobiastmp/C6280_S7_R1_001.fastq.gz -r2 /tmp/tobiastmp/C6280_S7_R2_001.fastq.gz -o /tmp/tobiastmp/outtest5 -ul 3 -sl 2 -dual
cp -r /tmp/tobiastmp/outtest5 .

