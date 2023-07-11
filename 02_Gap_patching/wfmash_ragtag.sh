#!/bin/sh

prefix=$1
ref=$2
region=$3
mkdir "$region"
cd "$region"
wfmash "$ref" "$prefix".fa > "$region".paf 
mkdir ragtag_output
cd ragtag_output
ln -s ../"$region".paf ragtag.patch.asm.paf
cd ..
ln -s "$ref" ref.fasta
ln -s "$prefix".fa query.fasta

ragtag.py patch -f 10000 --remove-small ref.fasta query.fasta


