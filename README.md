# SPART
SPART, a Semi-automated pipeline for assembling reference sequence of telomere-to-telomere (T2T). 
<img src="https://github.com/liushoucheng/SPART/blob/main/pic/pipeline.jpg" width="60%">

**See [tutorial]( https://spart1.readthedocs.io/en/latest/) for more details.** 
## Table of Contents

- [Quick install and start](#started)
  - [Install](#Install)
  - [Dependencies](#Dependencies)
  - [Running pipeline with snakemake](#pipe)
  - [Output files](#Output)
- [Run step by step](#step)
  - [00_Contig screen](#00_Contig)
  - [01_Contig scaffolding](#01_Contig)
  - [02_Gap patching](#02_Gap)
  - [03_Polishing](#03_Polishing)
  - [04_Evaluation](#04_Evaluation)
  - [05_Annotation](#05_Annotation)

## <a name="started"></a>Quick install and start
### <a name="Install"></a>Install
```sh
git clone https://github.com/liushoucheng/SPART.git
cd SPART
conda env create -f SPART.yaml
conda activate spart
```
### <a name="Dependencies"></a>Dependencies

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

### <a name="pipe"></a>Using snakemake to run the pipeline can be assembled to the chromosome level but may contain gaps that require the rest to be done manually.(Exclude Verkko,Bionano DLS Map,Telomere determination and patch,Centromeric region analysis,Variant calls and Evaluation):
* [Download the example in SPART/example/]( https://gofile.me/77wE8/Vj6Vlp1LK)
* [Download the digest_genome.py of HiC-Pro in SPART/]( https://github.com/nservant/HiC-Pro/blob/master/bin/utils/digest_genome.py)
```sh
# Replace SPART_PATH with the current working directory
sed -i "s#^ SPART_PATH# ${PWD}#g" conf_ck.yaml
# HiC enzyme
HiC_enzyme=" GATC"
# Replace hic_sca_enzyme with the value stored in the HiC_enzyme variable
sed -i "s#^ hic_sca_enzyme# ${HiC_enzyme}#g" conf_ck.yaml
# Ligation site sequence used for reads trimming. Depends on the fill in strategy. Example: AAGCTAGCTT
HiC_ligation_site=" GATCGATC"
sed -i "s#^ hic_sca_ligation_site# ${HiC_ligation_site}#g" conf_ck.yaml #Replace hic_sca_ligation_site with the value stored in the HiC_ligation_site variable
# This process uses the centos 7.6 operating system, slurm job scheduling system, please modify your SPART/clust.json according to the cluster situation.
# This process requires the use of HiC-Pro, please add it to the environment before running.
snakemake -s SPART.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs $threads --rerun-incomplete --restart-times 1 -np --rulegraph |dot -Tpng > rule.png #Running pipeline with snakemake
# configfile:The config file can be used to define a dictionary of configuration parameters and their values.
# cluster-config:A JSON or YAML file that defines the wildcards used in 'cluster'for specific rules.
```
<div align=center>
<img src="https://github.com/liushoucheng/SPART/blob/main/pic/rule.png" width="15%">
</div>

### <a name="Output"></a>Output files
please see the complete [documentation](https://github.com/liushoucheng/SPART/tree/main/example).

## <a name="step"></a>Run step by step

### <a name="00_Contig"></a>00_Contig screen
```sh
HiFi_reads=# file names of HiFi reads
ONT_reads=# file names of Ultra-Long reads
thread=# number of threads
memory=# Specify the upper limit on memory to use
output_prefix=# prefix of output files
mitochondrion=# mitochondrion fasta
chloroplast=# chloroplast fasta
ref=# Sequences of mitochondria and chloroplasts need to be removed
# Fastp :was used to filter adapter sequences, primers and other low quality sequence from raw sequencing reads.
SPART/00_Contig_screen/fastp.sh $HiFi_reads $ONT_reads
# Hifiasm
SPART/00_Contig_screen/hifiasm.sh $HiFi_reads $ONT_reads $output_prefix $thread
# Verkko
SPART/00_Contig_screen/verkko.sh $output_prefix $HiFi_reads $ONT_reads $threads $memory
# Flye
SPART/00_Contig_screen/flye.sh $ONT_reads $output_prefix $threads
# Remove mitochondrion && chloroplast
SPART/00_Contig_screen/rm_mt_cp.sh $mitochondrion $chloroplast $ref $threads
```
### <a name="01_Contig"></a>01_Contig scaffolding
```sh
threads=# Nominal threads per Node, without overloading (non-zero value will override -T -Tp -Te -TJ)
bnx=# Input molecule (.bnx) file, required
ref_cmap=# Reference file (must be .cmap), to compare resulting contigs
prefix=# Location of output files root directory, required, will be created if does not exist; if does exist, will overwrite contents
xml=# Read XML file for parameters
Bio_dir=# Location of executable files (RefAligner and Assembler, required)
cluster_xml=# Run on cluster, read XML file for submission arguments (optional--will not use cluster submission if absent)
ref=# Input NGS FASTA
bio_camp=# Input BioNano CMAP
merge_xml=# Merge configuration file
RefAligner=# RefAligner program
hicpro_data=# input data folder; Must contains a folder per sample with input files
hicpro_config=# configuration file for Hi-C processing
hicpro_outdir=# output folder
enzyme=# restriction enzyme cutting sites
#### Bionano
SPART/01_Contig_scaffolding/Bionano_DLS_map.sh $threads $bnx $ref_cmap $prefix $xml $Bio_dir $cluster_xml $ref $bio_camp $merge_xml $RefAligner
#### Hi-C
# hic-pro
SPART/01_Contig_scaffolding/HiC-Pro.sh $ref $prefix $hicpro_data $hicpro_config $hicpro_outdir
# yahs
SPART/01_Contig_scaffolding/yahs.sh $enzyme $ref $bed/bam/bin $profix
```
### <a name="02_Gap"></a>02_Gap patching
```sh
query=# query fasta file (uncompressed or bgzipped)
ref=# target fasta file (uncompressed or bgzipped)
region=# output directory
SPART/02_Gap_patching/wfmash_ragtag.sh $query $ref $region
```
#### Manual operation
```sh
cd ragtag_output
perl SPART/02_Gap_patching/paf_filter.pl -i ragtag.patch.debug.filtered.paf -minlen 10000000 -iden 0.5
```
**Manually editing the ragtag.patch.debug.filtered.paf file.Keep the high-quality contig and preserve the location of the only high confidence match in ragtag.patch.debug.filtered.paf that matches the sequence at both ends of the gap.**
```sh
perl SPART/02_Gap_patching/renameagp.pl -i ragtag.patch.ctg.agp -i1 ragtag.patch.debug.filtered.paf -start seq00000000 -end seq00000001 -o test.agp
```
**Test.agp is merged into ragtag.patch.agp and fasta is generated.**

#### e.g.
```sh
# make joins and fill gaps in target.fa using sequences from query.fa
cd SPART/example
ragtag.py patch -i 0.99 --remove-small -q 10 --debug -u --aligner minimap2 -t 128 --mm2-params "-x asm20 -I1G -t 128" reference1A.fasta query1A.fasta
# filter
cd ragtag_output
perl SPART/02_Gap_patching/paf_filter.pl -i ragtag.patch.debug.filtered.paf -minlen 10000000 -iden 0.5
# Manually editing the ragtag.patch.debug.filtered.paf_fiter.paf file.Keep the high-quality contig and preserve the location of the only high confidence match in ragtag.patch.debug.filtered.paf_fiter.paf that matches the sequence at both ends of the gap.
less ragtag.patch.debug.filtered.paf_fiter.paf
qseq00000000    600453479       27150   3999147 +       seq00000001     3972000 4       3971997 2266668 3972018 60
qseq00000000    600453479       4038251 35116708        +       seq00000002     597339226       17      31075089        17568679        31079144        60
# gain agp
perl SPART/02_Gap_patching/renameagp.pl -i ragtag.patch.ctg.agp -i1 ragtag.patch.debug.filtered.paf_fiter.paf -start seq00000001 -end seq00000002 -o test.agp
less -S ragtag.patch.agp
chr1A_RagTag_MOD_MOD	1	2046621	1	W	seq00000000	1	2046621	+
chr1A_RagTag_MOD_MOD	2046622	2046821	2	N	200	scaffold	yes	align_genus
chr1A_RagTag_MOD_MOD	2046822	6018821	3	W	seq00000001	1	3972000	+
chr1A_RagTag_MOD_MOD	6018822	6019021	4	N	200	scaffold	yes	align_genus
chr1A_RagTag_MOD_MOD	6019022	603358247	5	W	seq00000002	1	597339226	+
# Test.agp is merged into ragtag.patch.agp and fasta is generated.
less -S ragtag.patch.agp
scf00000000	1	2046621	1	W	seq00000000	1	2046621	+
scf00000000	2046622	2046821	2	N	200	scaffold	yes	align_genus
scf00000000	2046822	6018821	3	W	seq00000001	1	3972000	+
scf00000000	6018822	6057905	4	W	qseq00000000	3999151	4038234	+
scf00000000	6057906	603397131	5	W	seq00000002	1	597339226	+
ragtag_agp2fa.py ragtag.patch.agp ragtag.patch.comps.fasta > ragtag.patch.fasta
```
#### telomere patching
We used _submit_telomere.sh in ONT reads >100kb.ONT reads with telomere sequence mapping to this locus based on minimap2 alignments were manually identified. The longest was selected as template , all others aligned to it and polished with Medaka:
```sh
medaka -v -i ONT_tel_reads.fasta -d longest_ont_tel.fasta -o ont_tel_medaka.fasta
```
Telomere signal in all HiFi reads was identified with the commands:
```sh
_submit_telomere.sh hifi_reads.fasta
```
Additional HiFi reads were recruited from a manual analysis. We looked for trimmed tips that could extend. All reads had telomere signal and were aligned to the medaka consensus and polished with Racon with the commands:
```sh
minimap2 -t16 -ax map-pb ont_tel_medaka.fasta hifi_tel.fasta > medaka.sam
racon hifi_tel.fasta medaka.sam ont_tel_medaka.fasta > racon.fasta
```
Finally, the polished result was patched into the assembly with ragtag patch or manually patched.
##### Citation
https://github.com/marbl/CHM13-issues/blob/main/error_detection.md.
#### Centromeric region analysis
```sh
workdir=# work directory
FASTA=# target fasta file (uncompressed or bgzipped)
prefix=# prefix of output files
CHIP1_treatment=# Treatment (pull-down) file(s).
CHIP2_treatment=# Treatment (pull-down) file(s).
threads=# number of threads
CHIP1_control=# Control (input) file(s)
CHIP2_control=# Control (input) file(s)
SPART/02_Gap_patching/Centromeric_region_analysis.sh $workdir $FASTA $prefix $CHIP1_treatment $CHIP2_treatment $threads $CHIP1_control $CHIP2_control
```
### <a name="03_Polishing"></a>03_Polishing
```sh
# Use singularity and docker to download google_deepvariant_latest-gpu.sif and kishwars/pepper_deepvariant:r0.8-gpu respectively and modify the cluster-config and configfile in snakemake
workdir=# work directory
ref=# target fasta file (uncompressed or bgzipped)
threads=# number of threads
SPART/03_Polishing/calsv_snv.sh $workdir $ref $threads
```
### <a name="04_Evaluation"></a>04_Evaluation
```sh
ref=# target fasta file (uncompressed or bgzipped)
prefix=# prefix of output files
query=# query fasta file (uncompressed or bgzipped)
threads=# number of threads
partition=# your cluster partition
bac_reads=# bac reads
ref_chr=# target chromosome fasta file (uncompressed or bgzipped)
protein=# target protein fasta file
name=# output file name
gff3# target gff file
#### BUSCO
SPART/04_Evaluation/BUSCO.sh $ref $prefix
#### mapping rates & coverages
SPART/04_Evaluation/mapping_rates_coverages.sh hybrid_bam single_bam ont_bam
#### LTR
SPART/04_Evaluation/ltr.sh $ref $prefix
#### QV
SPART/04_Evaluation/qv.sh $query $ref
#### BACs
SPART/04_Evaluation/bac.sh $bac_reads $ref_chr
### Addition
SPART/04_Evaluation/while.sh $threads $partition $ref $query
### Analysis of synteny
SPART/04_Evaluation/synteny.sh $protein $name $gff3
```
### <a name="05_Annotation"></a>05_Annotation
#### RNA-seq
Detect adapter
```sh
fastp --detect_adapter_for_pe -w ${threads} -i ${RNAseq1} -I ${RNAseq2} -o ${RNAseq1_clean} -O ${RNAseq2_clean} --json ${output}.json --html ${output}.html
```
Build genome index
```sh
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${Output Dir} --genomeFastaFiles ${genome} --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile ${annotations} --limitGenomeGenerateRAM 40000000000 --sjdbOverhang 149 --sjdbFileChrStartEnd ${genomic coordinates} --limitSjdbInsertNsj 1854820
```
Mapping to genome
```sh
STAR --runThreadN ${threads} --genomeDir ${Output Dir} --readFilesIn ${RNAseq1_clean} ${RNAseq2_clean} --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile ${annotations} --outFileNamePrefix "$profix" --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterType BySJout --outSAMunmapped Within --outFilterMultimapNmax 20 --outSAMstrandField intronMotif --outFilterMismatchNoverLmax 0.02 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 10000 --alignMatesGapMax 100000 --sjdbScore 1 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --limitSjdbInsertNsj 1854820
```
Assembly and merge
```sh
stringtie -j 2 -c 2 -m 150 -f 0.3 -G ${reference annotation} -l rna-seq -t -p ${threads} -l "$profix" -A "$profix"gene_abund.tab -C "$profix"cov_refs.gtf -o "$profix".gtf "$profix"Aligned.sortedByCoord.out.bam
stringtie --merge -p 96 -m 150 -c 10 -G ${reference annotation} -l rna_merge -o rna_all.gtf { gtf_list | strg1.gtf ...}
```
#### ISO-seq
Build genome index
```sh
minimap2 -t 96 -I 16G -d $mmi $genome
```
Mapping to genome
```sh
flair 123 --mm2_args=-I15g,-axsplice:hq,-uf,-secondary=no -g $genome -r $iso_seq --mm_index $mmi -f $gtf -o flair.output --temp_dir temp_flair --stringent --no_gtf_end_adjustment --check_splice --generate_map --trust_end -t 96 --annotation_reliant generate --junction_bed $stringtie.bed
```
# Contacts

shoucheng Liu (liusc_work@163.com)

xiaopeng Li (xiaopeng.li@pku-iaas.edu.cn)
