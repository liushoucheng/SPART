# 00_Contig screen
## Fastp :was used to filter adapter sequences, primers and other low quality sequence from raw sequencing reads.
SPART/00_Contig_screen/fastp.sh $HiFi_reads $ONT_reads
## Hifiasm
SPART/00_Contig_screen/hifiasm.sh $HiFi_reads $ONT_reads $output_prefix
## Verkko
SPART/00_Contig_screen/verkko.sh $output_prefix $HiFi_reads $ONT_reads $threads $memory
## Flye
SPART/00_Contig_screen/flye.sh $ONT_reads $output_prefix $threads
## Remove MT & CP
SPART/00_Contig_screen/rm_mt_cp.sh $mitochondrion $chloroplast $ref
# 01_Contig scaffolding
## Bionano
SPART/01_Contig_scaffolding/Bionano_DLS_map.sh threads bnx ref_cmap prefix xml Bio_dir cluster_xml ref bio_camp merge_xml RefAligner
## Hi-C
SPART/01_Contig_scaffolding/HiC-Pro.sh ref ref_prefix hicpro_data hicpro_config hicpro_outdir

SPART/01_Contig_scaffolding/yahs.sh enzyme ref bed/bam/bin profix
# 02_Gap patching
SPART/02_Gap_patching/wfmash_ragtag.sh prefix ref region

## Manual operation

cd ragtag_output

perl SPART/02_Gap_patching/paf_filter.pl -i ragtag.patch.debug.filtered.paf -minlen 10000000 -iden 0.5

**Manually editing the ragtag.patch.debug.filtered.paf file.Keep the high-quality contig and preserve the location of the only high confidence match in ragtag.patch.debug.filtered.paf that matches the sequence at both ends of the gap.**

perl SPART/02_Gap_patching/renameagp.pl -i ragtag.patch.ctg.agp -i1 ragtag.patch.debug.filtered.paf -start seq00000000 -end seq00000001 -o test.agp

**Test.agp is merged into ragtag.patch.agp and fasta is generated.**

## telomere patching
We used _submit_telomere.sh in ONT reads >100kb.ONT reads with telomere sequence mapping to this locus based on minimap2 alignments were manually identified. The longest was selected as template , all others aligned to it and polished with Medaka:

medaka -v -i ONT_tel_reads.fasta -d longest_ont_tel.fasta -o ont_tel_medaka.fasta

Telomere signal in all HiFi reads was identified with the commands:

_submit_telomere.sh hifi_reads.fasta

Additional HiFi reads were recruited from a manual analysis. We looked for trimmed tips that could extend. All reads had telomere signal and were aligned to the medaka consensus and polished with Racon with the commands:

minimap2 -t16 -ax map-pb ont_tel_medaka.fasta hifi_tel.fasta > medaka.sam

racon hifi_tel.fasta medaka.sam ont_tel_medaka.fasta > racon.fasta

Finally, the polished result was patched into the assembly with ragtag patch or manually patched.
### Citation
https://github.com/marbl/CHM13-issues/blob/main/error_detection.md.
## Centromeric region analysis

SPART/02_Gap_patching/Centromeric_region_analysis.sh workdir FASTA INDEX prefix CHIP1 CHIP2 threads

# 03_Polishing
SPART/03_Polishing/calsv_snv.sh workdir ref threads
# 04_Evaluation
## BUSCO
SPART/04_Evaluation/BUSCO.sh ref prefix
## mapping rates & coverages
SPART/04_Evaluation/mapping_rates_coverages.sh hybrid_bam single_bam ont_bam
## LTR
SPART/04_Evaluation/ltr.sh ref prefix
## QV
SPART/04_Evaluation/qv.sh query ref
## BACs
SPART/04_Evaluation/bac.sh bac_reads ref_chr
