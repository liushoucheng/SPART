#!/bin/sh

mkdir blastdb
mkdir blastresult
protein=$1
name=$2
gff3=$3

makeblastdb -in $protein -dbtype prot -out ./blastdb/${name}
blastp -query $protein -db ./blastdb/${name} -out ./blastresult/${name}.blast -num_threads 52 -outfmt 6 -evalue 1e-10 -num_alignments 5

awk -vFS="\t" -vOFS="\t" '{if($3=="mRNA"){match($9,/ID=([^;]+)/,a);sub(/ID=/,"",a[0]);print $1,a[0],$4,$5}}' ${gff3} > ./blastresult/${name}.gff
cd blastresult
MCScanX ./${name}
