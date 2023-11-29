#!/bin/sh

workdir=$1
FASTA=$2
prefix=$3
CHIP1=$4
CHIP2=$5
threads=$6
CHIP1c=$7
CHIP2c=$8

#treatment
SAM="${workdir}/ref_chip.sam"
BAM="${workdir}/ref_chip.bam"
BAMSORT="${workdir}/ref_chip_sort.bam"
BAMFILTER="${workdir}/ref_chip_sort_filter.bam"
CHIP1TRIM="${workdir}/CHIP_1.trim.fastq"
CHIP2TRIM="${workdir}/CHIP_2.trim.fastq"
BINS="${workdir}/"$prefix".genome.size.100kb"

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
samtools index -c -@ $threads ${BAMFILTER}

#control
SAMc="${workdir}/ref_control.sam"
BAMc="${workdir}/ref_control.bam"
BAMSORTc="${workdir}/ref_control_sort.bam"
BAMFILTERc="${workdir}/ref_control_sort_filter.bam"
CHIP1TRIMc="${workdir}/ref_control.trim.fastq"
CHIP2TRIMc="${workdir}/ref_control.trim.fastq"
INTERSECTc="${workdir}/"$profix".intersect.bed"

# adapter
#fastp --in1 ${CHIP1c} --in2 ${CHIP2c}  --out1 ${CHIP1TRIMc} --out2 ${CHIP2TRIMc} --thread 52 

## RUN HISAT2
hisat2 -p $threads -x "$profix" -1 ${CHIP1TRIMc} -2 ${CHIP2TRIMc} -S ${SAMc}

## convert to BAM
samtools view -b -S -@ $threads -o ${BAMc} ${SAMc}

## sort
samtools sort -m 4G -@ $threads -o ${BAMSORTc} ${BAMc}

## filter
samtools view -@ $threads -q 30 -o ${BAMFILTERc} ${BAMSORTc}
samtools index -c -@ $threads ${BAMFILTERc}

epic2 -t ${BAMFILTER} -c ${BAMFILTERc} --chromsizes  -o "CENH3.bed" --bin-size 25000 --mapq 30 --gaps-allowed 4
