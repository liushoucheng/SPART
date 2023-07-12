# SPART
SPART, a Semi-automated pipeline for assembling reference sequence of telomere-to-telomere (T2T). 
<img width="703" alt="image" src="https://github.com/liushoucheng/SPART/assets/50602960/254b12f0-f3c7-4201-b9d2-f4a49876dd66">

## Quick install and start
### Install
git clone https://github.com/liushoucheng/SPART.git

cd SPART

conda env create -f SPART.yaml

cd ..

### Dependencies

Bionano DLS map https://bionano.com

HiC-Pro v3.1.0 https://github.com/nservant/HiC-Pro

DeepVariant v1.5.0-gpu https://github.com/google/deepvariant

PEPPER-Margin-DeepVariant v0.8-gpu https://github.com/kishwarshafin/pepper

hap.py v0.3.15 https://github.com/Illumina/hap.py

vcf_merge_t2t.py https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/polishing_merge_script/vcf_merge_t2t.py

### Running pipeline with snakemake:

snakemake -s SPART.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs $threads --rerun-incomplete --restart-times 1

configfile:The config file can be used to define a dictionary of configuration parameters and their values.

cluster-config:A JSON or YAML file that defines the wildcards used in 'cluster'for specific rules.

## Run step by step:

## 00_Contig screen
### Fastp :was used to filter adapter sequences, primers and other low quality sequence from raw sequencing reads.
SPART/00_Contig_screen/fastp.sh $HiFi_reads $ONT_reads
### Hifiasm
SPART/00_Contig_screen/hifiasm.sh $HiFi_reads $ONT_reads $output_prefix
### Verkko
SPART/00_Contig_screen/verkko.sh $output_prefix $HiFi_reads $ONT_reads $threads $memory
### Flye
SPART/00_Contig_screen/flye.sh $ONT_reads $output_prefix $threads
### Remove MT & CP
SPART/00_Contig_screen/rm_mt_cp.sh $mitochondrion $chloroplast $ref
## 01_Contig scaffolding
### Bionano
SPART/01_Contig_scaffolding/Bionano_DLS_map.sh threads bnx ref_cmap prefix xml Bio_dir cluster_xml ref bio_camp merge_xml RefAligner
### Hi-C
SPART/01_Contig_scaffolding/HiC-Pro.sh ref ref_prefix hicpro_data hicpro_config hicpro_outdir

SPART/01_Contig_scaffolding/yahs.sh enzyme ref bed/bam/bin profix
## 02_Gap patching
SPART/02_Gap_patching/wfmash_ragtag.sh prefix ref region
## 03_Polishing
SPART/03_Polishing/calsv_snv.sh workdir ref threads
## 04_Evaluation
### BUSCO
SPART/04_Evaluation/BUSCO.sh ref prefix
### mapping rates & coverages
SPART/04_Evaluation/mapping_rates_coverages.sh hybrid_bam single_bam ont_bam
### LTR
SPART/04_Evaluation/ltr.sh ref prefix
### QV
SPART/04_Evaluation/qv.sh query ref
### BACs
SPART/04_Evaluation/bac.sh bac_reads ref_chr
