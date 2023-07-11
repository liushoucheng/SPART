#!/bin/sh

ont=$1
hifi=$2

fastp -w 16 -i $hifi -o hifi_clean_data.fastq
fastp -q 10 -l 100000 -w 16 -i $ont -o ont_clean_data.fastq
