#!/bin/sh


hybrid_bam=$1
single_bam=$2
ont_bam=$3


samtools flagstat -@ 128 $hybrid_bam > hybrid_bam.flagstat
samtools coverage -o hybrid_bam.cov hybrid_bam
samtools flagstat -@ 128 $single_bam > single_bam.flagstat
samtools coverage -o single_bam.cov single_bam
samtools flagstat -@ 128 $ont_bam > ont_bam.flagstat
samtools coverage -o ont_bam.cov ont_bam
