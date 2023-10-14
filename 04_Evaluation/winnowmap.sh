#!/bin/sh

chr=$1
mkdir "$chr"
cd "$chr"
ref="$2".split/Chr"$chr".fa
query="$3".split/Chr"$chr".fa
samtools faidx $ref
cut -f 1,2 "$ref".fai > Chr"$chr".sizes

meryl count k=27 output merylDB_"$chr" ${ref}
meryl print greater-than distinct=0.9998 merylDB_"$chr" > repetitive_"$chr".txt
split_fa ${query} > split.fa

winnowmap -W repetitive_"$chr".txt -ax asm20 -K 1500M -k 27 -w 18 -t 52 -H --MD ${ref} split.fa > chr"$chr".sam

k8 paftools.js sam2paf -p chr"$chr".sam > chr"$chr".paf
cat chr"$chr".paf |awk '{if ($12 > 0) print $6"\t"$8"\t"$9}' |bedtools sort -i - |bedtools merge -i - |bedtools complement -i - -g Chr"$chr".sizes > Chr"$chr".bed

