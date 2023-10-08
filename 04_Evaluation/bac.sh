#!/bin/sh

item=$1
chr=$2
blastn -query $item -db ./"$chr"index -evalue 1e-6 -outfmt "6 qseqid qlen sseqid qstart qend sstart send pident length nident qcovs" -num_threads 128 -out "$1".bed
#-outfmt "6 qseqid qlen sseqid qstart qend sstart send pident length nident qcovs" 
sed -n '1p' "$1".bed > "$1".txt
mv "$1".txt ./
cat *.fasta.txt > all_bac.bed
#for chr in Chr3A Chr3B Chr5A Chr5D; do makeblastdb -in part_"$chr".fasta -dbtype nucl -parse_seqids -out ./"$chr"index; for item in dir bac.fasta.split/"$chr"/*.fasta; do sbatch --job-name="$chr"blast --partition= --cpus-per-task=128 blast.sh $item $chr; done; done