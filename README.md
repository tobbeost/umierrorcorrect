# umierrorcorrect

Pipeline for analyzing  barcoded amplicon sequencing data with Unique molecular identifiers (UMI)

Installation
------------

To install the UMI-errorcorrect pipeline, open a terminal and type the following:

```bash
wget https://github.com/tobbeost/umierrorcorrect/archive/v0.1.tar.gz
pip install v0.1.tar.gz
```

After installation, try to run the pipeline:

```bash
run_umierrorcorrect.py -h
```

Dependencies
------------

Umi-errorcorrect runs using Python and  requires the following programs/libraries to be installed:

Python-libraries (should be installed automatically):

    pysam (v 0.8.4 or greater)

External programs:

    bwa (bwa mem command is used)
    pigz

Install the external programs and add them to the path.


Usage
-----

Example syntax for running the whole pipeline:

    run_umierrorcorrect.py -r1 read1.fastq.gz -r2 read2.fastq.gz -ul umi_length -sl spacer_length -r reference_fasta_file.fasta -o output_directory

The ``run_umierrorcorrect.py`` pipeline performs the following steps:

- Preprocessing of fastq files (remove the UMI and spacer and puts the UMI in the header)
- Mapping of preprocessed fastq reads to the reference genome
- Perform UMI clustering, then error correcion of each UMI cluster
- Create consensus reads (one representative read per UMI cluster written to a BAM file)
- Create a consensus output file (collapsed counts per position)

It is also to possible to run the pipeline step-by-step.

To see the options for each step, type the following:

```bash
preprocess.py -h
run_mapping.py -h
umi_error_correct.py -h
get_consensus_statistics.py -h
filter_bam.py -h
filter_cons.py -h
```


