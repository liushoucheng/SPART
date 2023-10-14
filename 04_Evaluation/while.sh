#!/bin/sh

threads=$1
partition=$2
ref=$3
query=$4
seqkit split -i $ref
seqkit split -i $query
for num in {1..7}
do
for chr in A B D
do
	sbatch --job-name="$num""$chr" --partition=cuPartition --cpus-per-task="$threads" winnowmap.sh "$num""$chr" $ref $query
done
done
