# SPART
SPART, a Semi-automated pipeline for assembling reference sequence of telomere-to-telomere (T2T). 
![image](https://github.com/liushoucheng/SPART/blob/main/pic/pipeline.jpg)

See [tutorial]( https://spart1.readthedocs.io/en/latest/) for more details. 
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

### <a name="pipe"></a>Running pipeline with snakemake(Exclude Verkko,Bionano DLS Map,Telomere determination and patch,Centromeric region analysis,Variant calls and Evaluation):
```sh
sed -i "s#^ SPART_PATH# ${PWD}#g" conf_ck.yaml #Replace SPART_PATH with the current working directory
HiC_enzyme=" GATC" #HiC enzyme
sed -i "s#^ hic_sca_enzyme# ${HiC_enzyme}#g" conf_ck.yaml #Replace hic_sca_enzyme with the value stored in the HiC_enzyme variable
HiC_ligation_site=" GATCGATC" #Ligation site sequence used for reads trimming. Depends on the fill in strategy. Example: AAGCTAGCTT
sed -i "s#^ hic_sca_ligation_site# ${HiC_ligation_site}#g" conf_ck.yaml #Replace hic_sca_ligation_site with the value stored in the HiC_ligation_site variable
snakemake -s SPART.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs $threads --rerun-incomplete --restart-times 1 -np --rulegraph |dot -Tpng > rule.png #Running pipeline with snakemake
#configfile:The config file can be used to define a dictionary of configuration parameters and their values.
#cluster-config:A JSON or YAML file that defines the wildcards used in 'cluster'for specific rules.
```
<div align=center>
<img src="https://github.com/liushoucheng/SPART/blob/main/pic/rule.png">
</div>

### <a name="Output"></a>Output files
please see the complete [documentation]( https://github.com/liushoucheng/SPART/tree/main/exmple).

## <a name="step"></a>Run step by step

### <a name="00_Contig"></a>00_Contig screen
```sh
#### Fastp :was used to filter adapter sequences, primers and other low quality sequence from raw sequencing reads.
SPART/00_Contig_screen/fastp.sh $HiFi_reads $ONT_reads
#### Hifiasm
SPART/00_Contig_screen/hifiasm.sh $HiFi_reads $ONT_reads $output_prefix
#### Verkko
SPART/00_Contig_screen/verkko.sh $output_prefix $HiFi_reads $ONT_reads $threads $memory
#### Flye
SPART/00_Contig_screen/flye.sh $ONT_reads $output_prefix $threads
#### Remove mitochondrion && chloroplast
SPART/00_Contig_screen/rm_mt_cp.sh $mitochondrion $chloroplast $ref
```
### <a name="01_Contig"></a>01_Contig scaffolding
```sh
#### Bionano
SPART/01_Contig_scaffolding/Bionano_DLS_map.sh $threads $bnx $ref_cmap $prefix $xml $Bio_dir $cluster_xml $ref $bio_camp $merge_xml $RefAligner
#### Hi-C
SPART/01_Contig_scaffolding/HiC-Pro.sh $ref $ref_prefix $hicpro_data $hicpro_config $hicpro_outdir #hic-pro
SPART/01_Contig_scaffolding/yahs.sh $enzyme $ref $bed/bam/bin $profix #yahs
```
### <a name="02_Gap"></a>02_Gap patching
```sh
SPART/02_Gap_patching/wfmash_ragtag.sh $prefix $ref $region
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
SPART/02_Gap_patching/Centromeric_region_analysis.sh $workdir $FASTA $INDEX $prefix $CHIP1_treatment $CHIP2_treatment $threads $CHIP1_control $CHIP2_control
```
### <a name="03_Polishing"></a>03_Polishing
```sh
SPART/03_Polishing/calsv_snv.sh $workdir $ref $threads
```
### <a name="04_Evaluation"></a>04_Evaluation
```sh
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
