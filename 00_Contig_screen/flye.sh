#!/bin/sh

ont=$1
outdir=$2
threads=$3


flye --nano-hq $ont --read-error 0.1 -g 5.4g --asm-coverage 80 --scaffold --out-dir $outdir --threads $threads --no-alt-contigs