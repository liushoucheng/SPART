# SPART
SPART, a Semi-automated pipeline for assembling reference sequence of telomere-to-telomere (T2T). 
<img width="703" alt="image" src="https://github.com/liushoucheng/SPART/assets/50602960/254b12f0-f3c7-4201-b9d2-f4a49876dd66">

## dependence 
snakemake v7.21.0 https://snakemake.github.io
## 00_Contig screen
### fastp
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
