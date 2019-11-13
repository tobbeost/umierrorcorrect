# Docker documentation

If you have Docker installed, pull the Docker image from Docker hub.

```bash
docker pull tobiasosterlund/umierrorcorrect
```

Download a reference genome fasta file and mount the reference directory and data directory (including fastq files and BED files) to the docker container:

```bash
docker run -v /path_to_reference_fasta_directory/:/references/ -v /path_to_data_directory/:/data/ -it tobiasosterlund/umierrorcorrect
```
To try to run the pipeline:

```bash
run_umierrorcorrect.py -h
```

Since the umierrorcorrect pipeline is using `bwa` for mapping of reads, a bwa-indexed reference genome is needed. Index the reference genome with the command `bwa index -a bwtsw /references/reference.fa`.

Change the path to the data directory and run the pipeline. 
Example syntax for running the whole pipeline:

```bash
cd /data
run_umierrorcorrect.py -r1 read1.fastq.gz -r2 read2.fastq.gz -ul umi_length -sl spacer_length -r /references/reference.fa -o output_directory
```

