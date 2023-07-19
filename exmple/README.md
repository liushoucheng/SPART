# Running pipeline with snakemake(Exclude Verkko,Bionano DLS Map,Telomere determination and patch,Centromeric region analysis,Variant calls and Evaluation):
sed -i "s#^ SPART_PATH# ${PWD}#g" conf_ck.yaml

snakemake -s SPART.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs $threads --rerun-incomplete --restart-times 1 -np --rulegraph |dot -Tsvg > rule.svg

configfile:The config file can be used to define a dictionary of configuration parameters and their values.

cluster-config:A JSON or YAML file that defines the wildcards used in 'cluster'for specific rules.
# Output files
**Workdir/fastp**:Filtered adapter sequences, primers and other low quality 
sequence from raw HiFi and ONT sequencing reads.

**Workdir/hifiasm_hybrid/hybrid.all.asm.p_ctg.fa**:Hifiasm generated a preliminary contig genome assembly.

**Workdir/flye/assembly.fasta**:Flye generated the ONT UL reads assembly.

**Workdir/hifiasm_hybrid/hybrid.remove_cp_mt.fa**:Contigs with at least 50% of their bases covered by alignments were considered to be chloroplast or mitochondria genome sequences and were removed from the assembly.

**Workdir/hic_hybrid/hic_hybrid.bam**:Hi-C data were classified as valid or invalid interaction pairs.

**Workdir/yahs_hybrid/yahs_hybrid.fa**:Only valid interaction pairs were retained for subsequent assembly and scaffolding into chromosomes.

**Workdir/patch_flye/patch_single_hybrid_flye.fa**:Assembly gaps in chromosome scaffolds were directly filled by the corresponding Flye.

**Workdir/patch_verkko/patch_single_hybrid_flye_verkko.fa**:Assembly gaps in chromosome scaffolds were directly filled by the corresponding Verkko.

**Workdir/hybrid/hybrid.bam**:Alignment data file between patch_single_hybrid_flye_verkko.fa and HiFi reads.

**Workdir/hybrid_hifi_pcr/pcr.bam**:Alignment data file between patch_single_hybrid_flye_verkko.fa and PCR-FREE reads.

**Workdir/hybrid_hifi_pcr/hybrid.bam**:Merged Workdir/hybrid/hybrid.bam and Workdir/hybrid_hifi_pcr/pcr.bam

**Workdir/ont_merge/q10l120k.bam**:Alignment data file between patch_single_hybrid_flye_verkko.fa and ONT reads.