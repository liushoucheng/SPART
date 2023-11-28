import os
import re
import sys
b={}
hifi_single={}
hifi_mix={}
e={}
d={}
HiFi_hybrid_all=config["HiFi_reads_merge"]
ONT_all=config["ONT_reads_merge"]
mitochondrion=config["mitochondrion"]
chloroplast=config["chloroplast"]
hic_hybrid_dir=config["hic_dir"]
SPART_dir=config["SPART_dir"]
hic_hybrid_enzyme=config["hic_enzyme"]
hic_enzyme_ligation_site=config["hic_enzyme_ligation_site"]
verkko_fa=config["verkko_assemble"]
pcrfree_hybrid_r1=config["pcrfree_r1"]
pcrfree_hybrid_r2=config["pcrfree_r2"]
google_deepvariant_latest_gpu_sif=config["google_deepvariant_latest-gpu_sif"]
W=config["WORKDIR"]
DIR=config["DIR"]
DIRont=config["DIRont"]
for dirs in os.listdir(DIR):
    b2 = dirs.split(".fastq")
    if ".fastq" in dirs:
        absPath = os.path.join(DIR, dirs)
        hifi_mix[b2[0]]=absPath

for dirs in os.listdir(DIRont):
    b2 = dirs.split(".fastq")
    if ".fastq" in dirs:
        absPath = os.path.join(DIRont, dirs)
        e[b2[0]]=absPath

rule final:
    input:
        W+"hybrid_hifi_pcr/hybrid.bam",
        W + "ont_merge/q10l120k.bam"

rule hifi_fastp:
    input:
        HiFi_hybrid_all
    output:
        W+"fastp/hybrid.fq"
    shell:
        "fastp -w 16 -i {input} -o {output}"

rule ont_fastp:
    input:
        ONT_all
    output:
        W+"fastp/ont.fq"
    shell:
        "fastp -q 10 -l 100000 -w 16 -i {input} -o {output}"

rule hifiasm:
    input:
        hifi=W+"fastp/hybrid.fq",
        ont=W+"fastp/ont.fq"
    output:
        W+"hifiasm_hybrid/hybrid.all.asm.p_ctg.fa"
    params:
        W+"hifiasm_hybrid"
    shell:
        """
        cd {params}
        hifiasm -o hybrid.all.asm --primary -t 96 --ul {input.ont} -k 63 {input.hifi}
        awk '/^S/{{print ">"$2;print $3}}' hybrid.all.asm.p_ctg.gfa > {output}
        """

rule flye:
    input:
        W+"fastp/ont.fq"
    output:
        W + "flye/assembly.fasta"
    params:
        W
    shell:
        """
        cd {params}
        flye --nano-hq {input} --read-error 0.1 -g 5.4g --asm-coverage 80 --scaffold --out-dir flye --threads 96 --no-alt-contigs
        """

rule rm_mt_cp:
    input:
        hybrid=W+"hifiasm_hybrid/hybrid.all.asm.p_ctg.fa",
        mt=mitochondrion,
        cp=chloroplast
    output:
        W+"hifiasm_hybrid/hybrid.remove_cp_mt.fa"
    params:
        dir=W+"hifiasm_hybrid",
        workdir=SPART_dir
    shell:
        """
        cd {params.dir}
        minimap2 -t 96 -x asm5 {input.mt} {input.hybrid}> mitochondrion.paf
        minimap2 -t 96 -x asm5 {input.cp} {input.hybrid}> chloroplast.paf
        python {params.workdir}/gemma_los.py mitochondrion.paf > mitochondrion.txt
        python {params.workdir}/gemma_los.py chloroplast.paf > chloroplast.txt
        seqkit grep -v -f chloroplast.txt {input.hybrid} > wheat_remove_cp.fa
        seqkit grep -v -f mitochondrion.txt wheat_remove_cp.fa > {output}
        """

rule hicpro:
    input:
        hic=hic_hybrid_dir,
        ref=W+"hifiasm_hybrid/hybrid.remove_cp_mt.fa"
    output:
        W+"hic_hybrid/hic_hybrid.bam"
    params:
        dir=W+"hic_hybrid",
        prefix="hybrid.remove_cp_mt",
        spart_dir=SPART_dir,
        enzyme=hic_hybrid_enzyme,
        LIGATION_SITE=hic_enzyme_ligation_site
    shell:
        """
        cd {params.dir}
        ln -s {input.ref} ./
        bowtie2-build --large-index --threads 96 {params.prefix}.fa {params.prefix}
        samtools faidx {params.prefix}.fa
        awk '{{print $1 "\t" $2}}' {params.prefix}.fa.fai > genome_sizes.bed
        python {params.spart_dir}/digest_genome.py -r ^{params.enzyme} -o enzyme.bed {params.prefix}.fa
        makeblastdb -in {params.prefix}.fa -dbtype nucl -parse_seqids -out {params.prefix}
        cp {params.spart_dir}/01_Contig_scaffolding/hicpro_config.txt ./
        sed -i 's#^N_CPU = #N_CPU = 128#g' hicpro_config.txt
        sed -i 's#^BOWTIE2_IDX_PATH = #BOWTIE2_IDX_PATH = {params.dir}#g' hicpro_config.txt
        sed -i 's#^REFERENCE_GENOME = #REFERENCE_GENOME = {params.prefix}#g' hicpro_config.txt
        sed -i 's#^GENOME_SIZE = #GENOME_SIZE = {params.dir}/genome_sizes.bed#g' hicpro_config.txt
        sed -i 's#^GENOME_FRAGMENT = #GENOME_FRAGMENT = {params.dir}/enzyme.bed#g' hicpro_config.txt
        HiC-Pro -i {input.hic} -c hicpro_config.txt -o {params.dir}/result
        cd result/bowtie_results/bwt2/sample
        samtools sort -m 1500M -n -@ 96 HiC_hybrid.remove_cp_mt.bwt2pairs.bam > {params.dir}/hic_hybrid.bam
        """

rule yahs:
    input:
        bam=W+"hic_hybrid/hic_hybrid.bam",
        ref=W+"hifiasm_hybrid/hybrid.remove_cp_mt.fa"
    output:
        W + "yahs_hybrid/yahs_hybrid.fa"
    params:
        dir = W + "yahs_hybrid",
        prefix = "hybrid_bam",
        enzyme = hic_hybrid_enzyme
    shell:
        """
        cd {params.dir}
	samtools faidx {input.ref}
	samtools sort -@ 128 -o hic_hybrid_sort.bam {input.bam}
	samtools index hic_hybrid_sort.bam
        yahs -e {params.enzyme} {input.ref} hic_hybrid_sort.bam -o {params.prefix}
        cp {params.dir}/{params.prefix}_scaffolds_final.fa {output}
        """

rule patch_flye:
    input:
        single_hybrid=W + "yahs_hybrid/yahs_hybrid.fa",
        flye=W + "flye/assembly.fasta"
    output:
        W + "patch_flye/patch_single_hybrid_flye.fa"
    params:
        dir = W + "patch_flye",
        prefix = "single_hybrid_flye"
    shell:
        """
        cd {params.dir}
        samtools faidx {input.single_hybrid}
        samtools faidx {input.flye}
        wfmash {input.single_hybrid} {input.flye} > {params.prefix}.paf 
        mkdir ragtag_output
        cd ragtag_output
        ln -s ../{params.prefix}.paf ragtag.patch.asm.paf
        cd ..
        ragtag.py patch -f 10000 --remove-small {input.single_hybrid} {input.flye}
        cp {params.dir}/ragtag_output/ragtag.patch.fasta {output}
        """

rule patch_verkko:
    input:
        single_hybrid_flye=W + "patch_flye/patch_single_hybrid_flye.fa",
        verkko=verkko_fa
    output:
        ref=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        txt = W + "repetitive_k27.txt"
    params:
        dir = W + "patch_verkko",
        prefix = "single_hybrid_flye_verkko",
    shell:
        """
        cd {params.dir}
        samtools faidx {input.single_hybrid_flye}
        samtools faidx {input.verkko}
        wfmash {input.single_hybrid_flye} {input.verkko} > {params.prefix}.paf 
        mkdir ragtag_output
        cd ragtag_output
        ln -s ../{params.prefix}.paf ragtag.patch.asm.paf
        cd ..
        ragtag.py patch -f 10000 --remove-small {input.single_hybrid_flye} {input.verkko}
        cp {params.dir}/ragtag_output/ragtag.patch.fasta {output.ref}
        bwa-mem2 index {output.ref}
        meryl count k=27 output merylDB {output.ref}
        meryl print greater-than distinct=0.9998 merylDB > {output.txt}
        """

rule winnowmap_hifi:
    input:
        fq=W+"fastp/hybrid.fq",
        ref=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        txt = W + "repetitive_k27.txt"
    output:
        sam=W+"hifi_mix_winnowmap/{hifi_mix}_q40l15k.sam"
    benchmark:
        W+"benchmarks/hifi_mix_winnowmap/{hifi_mix}.benchmark.txt"
    shell:
        """
        winnowmap --MD -W {input.txt} -ax map-pb -H -K 1500M -k 27 -w27 -t32 {input.ref} {input.fq} > {output.sam}
        """

rule winnowmap_hifi_sort:
    input:
        W+"hifi_mix_winnowmap/{hifi_mix}_q40l15k.sam"
    output:
        W+"hifi_mix_sort/{hifi_mix}_q40l15k.bam"
    params:
        W + "patch_verkko/patch_single_hybrid_flye_verkko.fa.fai"
    benchmark:
        W + "benchmarks/hifi_mix_sort/{hifi_mix}.benchmark.txt"
    shell:
        "samtools view -@32 -bt {params} {input}|samtools sort -@32 -m1500M -O bam -o {output} -"

rule winnowmap_hifi_sort_filter:
    input:
        W+"hifi_mix_sort/{hifi_mix}_q40l15k.bam"
    output:
        W+"hifi_mix_sort_filter/{hifi_mix}_q40l15k.bam"
    benchmark:
        W + "benchmarks/hifi_mix_sort_filter/{hifi_mix}.benchmark.txt"
    shell:
        "samtools view -@32 -F0x104 -hb {input} > {output}"

rule winnowmap_hifi_sort_filter_merge:
    input:
        expand(W+"hifi_mix_sort_filter/{hifi_mix}_q40l15k.bam",hifi_mix=hifi_mix)
    output:
        W+"hybrid/hybrid.bam"
    benchmark:
        W + "benchmarks/hybrid/hybrid.benchmark.txt"
    shell:
        "samtools merge -@ 128 -l 0 {output} {input}"

rule pcr_free:
    input:
        fa=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        r1=pcrfree_hybrid_r1,
        r2=pcrfree_hybrid_r2
    output:
        W+"hybrid_hifi_pcr/pcr.bam"
    shell:
       "bwa-mem2.avx512bw mem -t 96 {input.fa} {input.r1} {input.r2}|samtools view -@ 96 -b -|samtools sort -@ 96 -m 30G -o {output} -"

rule winnowmap_hifi_filter_pcr_merge:
    input:
        hifi=expand(W+"hifi_mix_sort_filter/{hifi_mix}_q40l15k.bam",hifi_mix=hifi_mix),
        pcr=W+"hybrid_hifi_pcr/pcr.bam"
    output:
        W+"hybrid_hifi_pcr/hybrid.bam"
    benchmark:
        W + "benchmarks/hybrid_pcr/hybrid.benchmark.txt"
    shell:
        "samtools merge -@ 128 -l 0 {output} {input.hifi} {input.pcr}"

rule winnowmap_ont:
    input:
        fq=W+"fastp/ont.fq",
        ref=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        txt=W+"repetitive_k27.txt"
    output:
        W+"ont_winnowmap/{e}/{e}_q10l120k.sam"
    benchmark:
        W+"benchmarks/ont_winnowmap/{e}.benchmark.txt"
    shell:
        "winnowmap --MD -W {input.txt} -ax map-ont -H -K 1500M -k 27 -w27 -t32 {input.ref} {input.fq} > {output}"

rule winnowmap_ont_sort:
    input:
        W+"ont_winnowmap/{e}/{e}_q10l120k.sam"
    output:
        W+"ont_sort/{e}/{e}_q10l120k.bam"
    params:
        W + "patch_verkko/patch_single_hybrid_flye_verkko.fa.fai"
    benchmark:
        W + "benchmarks/ont_sort/{e}.benchmark.txt"
    shell:
        "samtools view -@32 -bt {params} {input}|samtools sort -@32 -m1500M -O bam -o {output} -"

rule winnowmap_ont_sort_filter:
    input:
        W+"ont_sort/{e}/{e}_q10l120k.bam"
    output:
        W+"ont_filter/{e}_q10l120k.bam"
    benchmark:
        W + "benchmarks/ont_filter/{e}.benchmark.txt"
    shell:
        "samtools view -@ 128 -F0x104 -hb {input} > {output}"

rule winnowmap_ont_sort_filter_merge:
    input:
        expand(W+"ont_filter/{e}_q10l120k.bam",e=e)
    output:
        W + "ont_merge/q10l120k.bam"
    benchmark:
        W + "benchmarks/ont_merge/benchmark.txt"
    shell:
        "samtools merge -@ 128 -o {output} {input}"
