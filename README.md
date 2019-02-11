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

Umi-errorcorrect requires the following programs/libraries to be installed:

Python-libraries (should be installed automatically):

    pysam (v 0.8.4 or greater)

External programs:

    bwa (bwa mem command is used)
    samtools
    pigz

Install the external programs and add them to the path.


