#!/bin/sh

snakemake -s bwa_winnowmap.py --cluster-config clust_align.json --configfile conf_ck_align.yaml --cluster '{cluster.account}' --jobs 40 --rerun-incomplete --restart-times 1
workdir=$1
input="$workdir"/input
output="$workdir"/output
mkdir -p $input
mkdir -p $output
mkdir -p "$output"/intermediate_results_dir/work
mkdir -p "$output"/intermediate_results_dir/temp
mkdir -p "$output"/pepper_deepvariant_output
ref=$2

a=${ref%.*}
b=${a##*/}
$(basename $ref .fa)
cp $ref $input
threads=$3
ref="$input"/"$b".fasta
mv ont_merge/q10l120k.bam $input
mv ont_merge/q10l120k.bam.csi $input
mv single_hifi_pcr/hybrid.bam $input
mv single_hifi_pcr/hybrid.bam.csi $input
mv hybrid_hifi_pcr/hybrid.bam $input
mv hybrid_hifi_pcr/hybrid.bam.csi $input

####dv

singularity run --cpus $threads --nv -B /home:/home -B /data:/data -B "$input":"$input" -B "$output":"$output" --workdir "$output"/intermediate_results_dir_hybrid/temp google_deepvariant_latest-gpu.sif /opt/deepvariant/bin/run_deepvariant --model_type "HYBRID_PACBIO_ILLUMINA" --ref "$ref" --reads "$input"/hybrid.bam --output_vcf "$output"/hybrid.vcf --num_shards $threads --intermediate_results_dir "$output"/intermediate_results_dir_hybrid

singularity run --cpus $threads --nv -B /home:/home -B /data:/data -B "$input":"$input" -B "$output":"$output" --workdir "$output"/intermediate_results_dir_single/temp google_deepvariant_latest-gpu.sif /opt/deepvariant/bin/run_deepvariant --model_type "HYBRID_PACBIO_ILLUMINA" --ref "$ref" --reads "$input"/single.bam --output_vcf "$output"/single.vcf --num_shards $threads --intermediate_results_dir "$output"/intermediate_results_dir_single

####pepper dv

docker run --ipc=host --gpus all -v /home:/home -v /data:/data -v "$input":"$input" -v "$output":"$output" kishwars/pepper_deepvariant:r0.8-gpu run_pepper_margin_deepvariant call_variant -b "$input"/q10l120k.bam -f $ref -o "$output"/pepper_deepvariant_output -g -p pep_dv_ont -t $threads --ont_r9_guppy5_sup

snakemake -s callsv_snv.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs 128 --rerun-incomplete --restart-times 1
