#!/bin/sh

hifi_reads=$1
ont_reads=$2
pre=$3

hifiasm -k 63 -o "$pre".asm -t 128 $hifi_reads --ul $ont_reads
