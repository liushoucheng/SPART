#!/bin/sh

threads=$1
bnx=$2
ref_cmap=$3
prefix=$4
xml=$5
Bio_dir=$6
cluster_xml=$7
ref=$8
bio_camp=$9
merge_xml=$10
RefAligner=$11

python pipelineCL.py -Tn $threads -i 5 -b $bnx -r $ref_cmap -l $prefix -e w -a $xml -t $Bio_dir -y -z --species-reference other -C $cluster_xml -F 1
perl hybridScaffold.pl -n $ref -b $bio_camp -u CTTAAG -c $merge_xml -r $RefAligner -o bio_hybrid -B 2 -N 2 -f