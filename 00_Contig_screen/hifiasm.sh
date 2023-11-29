#!/bin/sh

hifi_reads=$1
ont_reads=$2
pre=$3
thread=$4

hifiasm -k 63 -o "$pre".asm -t $thread $hifi_reads --ul $ont_reads
