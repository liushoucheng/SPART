#!/bin/sh

ref=$1
ref_prefix=$2
hicpro_data=$3
hicpro_config=$4
hicpro_outdir=$5

ln -s $ref ./"$ref_prefix".fa
bowtie2-build --large-index --threads 128 "$ref_prefix".fa "$ref_prefix"

samtools faidx "$ref_prefix".fa
awk '{print $1 "\t" $2}' "$ref_prefix".fa.fai > genome_sizes.bed

python ./HiC-Pro/bin/utils/digest_genome.py -r ^GATC -o wheat_DpnII.bed "$ref_prefix".fa
makeblastdb -in "$ref_prefix".fa -dbtype nucl -parse_seqids -out "$ref_prefix"

HiC-Pro -i $hicpro_data -c $hicpro_config -o $hicpro_outdir -p
