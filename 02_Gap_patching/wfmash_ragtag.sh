#!/bin/sh

query=$1
ref=$2
region=$3
mkdir "$region"
cd "$region"
wfmash "$ref" "$query" > "$region".paf 
mkdir ragtag_output
cd ragtag_output
ln -s ../"$region".paf ragtag.patch.asm.paf
cd ..
ln -s "$ref" ref.fasta
ln -s "$query" query.fasta

ragtag.py patch -f 10000 --remove-small ref.fasta query.fasta

_submit_telomere.sh ragtag_output/ragtag.patch.fasta #(https://github.com/VGP/vgp-assembly)
