#!/bin/sh

workdir=$1
FASTA=$2
INDEX=$3
prefix=$4
CHIP1=$5
CHIP2=$6
threads=$7
SAM="${workdir}/ref_chip.sam"
BAM="${workdir}/ref_chip.bam"
BAMSORT="${workdir}/ref_chip_sort.bam"
BAMFILTER="${workdir}/ref_chip_sort_filter.bam"
CHIP1TRIM="${workdir}/CHIP_1.trim.fastq"
CHIP2TRIM="${workdir}/CHIP_2.trim.fastq"
BINS="${workdir}/"$prefix".genome.size.100kb"
INTERSECT="${workdir}/"$prefix".intersect.bed"

faidx $FASTA -i chromsizes | bedtools makewindows -g - -w 100000 | awk -vFS="\t" -vOFS="\t" '{print $1,$2,$3}' | bedtools sort -i - > "$prefix".genome.size.100kb

ln -s ${FASTA} ./"$prefix".fasta
hisat2-build --large-index -a -p $threads ${workdir}/"$prefix".fasta ${workdir}/"$prefix"

# adapter
fastp --in1 ${CHIP1} --in2 ${CHIP2}  --out1 ${CHIP1TRIM} --out2 ${CHIP2TRIM} --thread $threads 

## RUN HISAT2
hisat2 -p $threads -x "$prefix" -1 ${CHIP1TRIM} -2 ${CHIP2TRIM} -S ${SAM} 

## convert to BAM
samtools view -b -S -@ $threads -o ${BAM} ${SAM}

## sort
samtools sort -m 4G -@ $threads -o ${BAMSORT} ${BAM}

## filter
samtools view -@ $threads -q 30 -o ${BAMFILTER} ${BAMSORT}

## intersect bed
bedtools intersect -a ${BINS} -b ${BAMFILTER} -c -sorted > ${INTERSECT}

python SPART/02_Gap_patching/chip-seq.py ${INTERSECT} > cen.bed
