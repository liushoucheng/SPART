# SPART
SPART, a Semi-automated pipeline for assembling reference sequence of telomere-to-telomere (T2T). 
![image](https://github.com/liushoucheng/SPART/blob/main/pic/pipeline.jpg)

## Quick install and start
### Install
git clone https://github.com/liushoucheng/SPART.git

cd SPART

conda env create -f SPART.yaml

conda activate spart

### Dependencies

List of tools assumed loadable or accessible with no path are:

* [Bionano DLS map]( https://bionano.com)

* [HiC-Pro v3.1.0]( https://github.com/nservant/HiC-Pro)

* [_submit_telomere.sh]( https://github.com/VGP/vgp-assembly/blob/master/pipeline/telomere/_submit_telomere.sh)

* [Medaka]( https://anaconda.org/bioconda/medaka)

* [racon]( https://anaconda.org/bioconda/racon)

* [hisat2]( https://github.com/DaehwanKimLab/hisat2)

* [DeepVariant v1.5.0-gpu]( https://github.com/google/deepvariant)

* [PEPPER-Margin-DeepVariant v0.8-gpu]( https://github.com/kishwarshafin/pepper)

* [hap.py v0.3.15]( https://github.com/Illumina/hap.py)

* [vcf_merge_t2t.py](https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/polishing_merge_script/vcf_merge_t2t.py)

### Running pipeline with snakemake(Exclude Verkko,Bionano DLS Map,Telomere determination and patch,Centromeric region analysis,Variant calls and Evaluation):

sed -i "s#^ SPART_PATH# ${PWD}#g" conf_ck.yaml

HiC_enzyme=" GATC"

sed -i "s#^ hic_sca_enzyme# ${HiC_enzyme}#g" conf_ck.yaml

HiC_ligation_site=" GATCGATC"

sed -i "s#^ hic_sca_ligation_site# ${HiC_ligation_site}#g" conf_ck.yaml

snakemake -s SPART.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs $threads --rerun-incomplete --restart-times 1 -np --rulegraph |dot -Tsvg > rule.svg

configfile:The config file can be used to define a dictionary of configuration parameters and their values.

cluster-config:A JSON or YAML file that defines the wildcards used in 'cluster'for specific rules.
![image](https://github.com/liushoucheng/SPART/blob/main/pic/rule.png)
#### Output files
**Workdir/fastp:**Filtered adapter sequences, primers and other low quality 
sequence from raw HiFi and ONT sequencing reads.

**Workdir/hifiasm_hybrid/hybrid.all.asm.p_ctg.fa:**Hifiasm generated a preliminary contig genome assembly.

**Workdir/flye/assembly.fasta:**Flye generated the ONT UL reads assembly.

**Workdir/hifiasm_hybrid/hybrid.remove_cp_mt.fa:**Contigs with at least 50% of their bases covered by alignments were considered to be chloroplast or mitochondria genome sequences and were removed from the assembly.

**Workdir/hic_hybrid/hic_hybrid.bam:**Hi-C data were classified as valid or invalid interaction pairs.

**Workdir/yahs_hybrid/yahs_hybrid.fa:**Only valid interaction pairs were retained for subsequent assembly and scaffolding into chromosomes.

**Workdir/patch_flye/patch_single_hybrid_flye.fa:**Assembly gaps in chromosome scaffolds were directly filled by the corresponding Flye.

**Workdir/patch_verkko/patch_single_hybrid_flye_verkko.fa:**Assembly gaps in chromosome scaffolds were directly filled by the corresponding Verkko.

**Workdir/hybrid/hybrid.bam:**Alignment data file between patch_single_hybrid_flye_verkko.fa and HiFi reads.

**Workdir/hybrid_hifi_pcr/pcr.bam:**Alignment data file between patch_single_hybrid_flye_verkko.fa and PCR-FREE reads.

**Workdir/hybrid_hifi_pcr/hybrid.bam:**Merged Workdir/hybrid/hybrid.bam and Workdir/hybrid_hifi_pcr/pcr.bam

**Workdir/ont_merge/q10l120k.bam:**Alignment data file between patch_single_hybrid_flye_verkko.fa and ONT reads.

## Run step by step:

### 00_Contig screen
#### Fastp :was used to filter adapter sequences, primers and other low quality sequence from raw sequencing reads.
SPART/00_Contig_screen/fastp.sh $HiFi_reads $ONT_reads
#### Hifiasm
SPART/00_Contig_screen/hifiasm.sh $HiFi_reads $ONT_reads $output_prefix
#### Verkko
SPART/00_Contig_screen/verkko.sh $output_prefix $HiFi_reads $ONT_reads $threads $memory
#### Flye
SPART/00_Contig_screen/flye.sh $ONT_reads $output_prefix $threads
#### Remove MT & CP
SPART/00_Contig_screen/rm_mt_cp.sh $mitochondrion $chloroplast $ref
### 01_Contig scaffolding
#### Bionano
SPART/01_Contig_scaffolding/Bionano_DLS_map.sh threads bnx ref_cmap prefix xml Bio_dir cluster_xml ref bio_camp merge_xml RefAligner
#### Hi-C
SPART/01_Contig_scaffolding/HiC-Pro.sh ref ref_prefix hicpro_data hicpro_config hicpro_outdir

SPART/01_Contig_scaffolding/yahs.sh enzyme ref bed/bam/bin profix
### 02_Gap patching
SPART/02_Gap_patching/wfmash_ragtag.sh prefix ref region
#### telomere patching
We used _submit_telomere.sh in ONT reads >100kb.ONT reads with telomere sequence mapping to this locus based on minimap2 alignments were manually identified. The longest was selected as template , all others aligned to it and polished with Medaka:

medaka -v -i ONT_tel_reads.fasta -d longest_ont_tel.fasta -o ont_tel_medaka.fasta

Telomere signal in all HiFi reads was identified with the commands:

_submit_telomere.sh hifi_reads.fasta

Additional HiFi reads were recruited from a manual analysis. We looked for trimmed tips that could extend. All reads had telomere signal and were aligned to the medaka consensus and polished with Racon with the commands:

minimap2 -t16 -ax map-pb ont_tel_medaka.fasta hifi_tel.fasta > medaka.sam

racon hifi_tel.fasta medaka.sam ont_tel_medaka.fasta > racon.fasta

Finally, the polished result was patched into the assembly with ragtag patch or manually patched.
##### Citation
https://github.com/marbl/CHM13-issues/blob/main/error_detection.md.
#### Centromeric region analysis

SPART/02_Gap_patching/Centromeric_region_analysis.sh workdir FASTA INDEX prefix CHIP1 CHIP2 threads

### 03_Polishing
SPART/03_Polishing/calsv_snv.sh workdir ref threads
### 04_Evaluation
#### BUSCO
SPART/04_Evaluation/BUSCO.sh ref prefix
#### mapping rates & coverages
SPART/04_Evaluation/mapping_rates_coverages.sh hybrid_bam single_bam ont_bam
#### LTR
SPART/04_Evaluation/ltr.sh ref prefix
#### QV
SPART/04_Evaluation/qv.sh query ref
#### BACs
SPART/04_Evaluation/bac.sh bac_reads ref_chr
