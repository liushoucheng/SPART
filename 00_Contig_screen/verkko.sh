#!/bin/sh


output=$1
HiFi=$2
ONT=$3
threads=$4
memory=$5

verkko -d $output --hifi $HiFi --nano $4 --threads $threads  --slurm --local-memory $memory --snakeopts "--max-jobs-per-second 10 --max-status-checks-per-second 0.5 --restart-times 1 --local-cores 128 --jobs 250" --base-k 1001 --window 971 --hifi-coverage 100 --slurm --sto-run 128 200 24 --mer-run 128 200 20000  --ovb-run 64 100 24 --ovs-run 16 35 24  --red-run 16 31 24 --mbg-run 32 0 20000 --utg-run 128 240 20000 --spl-run 128 240 20000 --ali-run 23 50 20000 --pop-run 128 240 20000 --utp-run 1 240 200000  --lay-run 1 240 20000  --sub-run 128 240 200000 --par-run 128 240 20000 --cns-run 24 0 20000
