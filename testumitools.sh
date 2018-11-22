#!/bin/bash
umi_tools extract --extract-method regex --bc-pattern="^(?P<umi_1>.{3})(?P<discard_1>.{2})" --bc-pattern2="^(?P<umi_2>.{3})(?P<discard_1>.{2})" --stdin /tmp/tobiastmp/C6280_S7_R1_001.fastq.gz --read2-in /tmp/tobiastmp/C6280_S7_R2_001.fastq.gz -S /tmp/tobiastmp/umi_1.fastq.gz --read2-out /tmp/tobiastmp/umi_2.fastq.gz
