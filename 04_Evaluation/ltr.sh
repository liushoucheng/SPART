#!/bin/sh

ref=$1
prefix=$2

LTR_FINDER_parallel -seq $ref -threads 96 -harvest_out -size 1000000

LTR_retriever -threads 96 -genome $ref -inharvest m2.1.7.fasta.finder.combine.scn -dnalib clariTeRep.fna -plantprolib protein.fasta
