#!/bin/sh

mt=$1
cp=$2
ref=$3
threads=$4
minimap2 -t $threads -x asm5 $mt $ref> mitochondrion.paf

minimap2 -t $threads -x asm5 $cp $ref> chloroplast.paf

python gemma_los.py mitochondrion.paf > mitochondrion.txt
python gemma_los.py chloroplast.paf > chloroplast.txt

seqkit grep -v -f chloroplast.txt $ref > wheat_remove_cp.fa

seqkit grep -v -f mitochondrion.txt wheat_remove_cp.fa > wheat_remove_cp_mt.fa