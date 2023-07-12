import os
import re
import sys
b={}
hifi_single={}
hifi_mix={}
e={}
d={}
HiFi_hybrid_all=sys.argv[1]
HiFi_single_all=sys.argv[2]
ONT_all=sys.argv[3]
mitochondrion=sys.argv[4]
chloroplast=sys.argv[5]
hic_hybrid_dir=sys.argv[6]
hic_single_dir=sys.argv[7]
SPART_dir=sys.argv[8]
hic_hybrid_enzyme=sys.argv[9]
hic_single_enzyme=sys.argv[10]
verkko_fa=sys.argv[11]
pcrfree_hybrid_r1=sys.argv[12]
pcrfree_hybrid_r2=sys.argv[13]
pcrfree_single_r1=sys.argv[14]
pcrfree_single_r2=sys.argv[15]
google_deepvariant_latest-gpu_sif=sys.argv[16]
W=config["WORKDIR"]
DIR=config["DIR"]
DIRs=config["DIRs"]
DIRont=config["DIRont"]
for dirs in os.listdir(DIR):
    b2 = dirs.split(".fastq")
    if ".fastq" in dirs:
        absPath = os.path.join(DIR, dirs)
        hifi_mix[b2[0]]=absPath
    elif ".bam" in dirs:
        a = 0

for x in range(1, 4):
    x=str(x)
    for dirs in os.listdir(DIRs+x):
        b2 = dirs.split(".fastq")
        if ".fastq" in dirs:
            absPath = os.path.join(DIRs+x, dirs)
            d[b2[0]]=absPath
        elif ".bam" in dirs:
            a = 0

for dirs in os.listdir(DIRont):
    b2 = dirs.split(".fastq")
    if ".fastq" in dirs:
        absPath = os.path.join(DIRont, dirs)
        e[b2[0]]=absPath
    elif ".bam" in dirs:
        a = 0

rule final:
    input:
        W+"single_hifi_pcr/hybrid.bam",
        W+"hybrid_hifi_pcr/hybrid.bam",
        W + "ont_merge/q10l120k.bam"
rule hybrid_fastp:
    input:
        HiFi_hybrid_all
    output:
        W+"fastp/hybrid.fq"
    shell:
        "fastp -w 16 -i {input} -o {output}"
rule single_fastp:
    input:
        HiFi_single_all
    output:
        W+"fastp/single.fq"
    shell:
        "fastp -w 16 -i {input} -o {output}"
rule ont_fastp:
    input:
        ONT_all
    output:
        W+"fastp/ont.fq"
    shell:
        "fastp -q 10 -l 100000 -w 16 -i {input} -o {output}"

rule hifiasm_hybrid:
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
        awk '/^S/{print ">"$2;print $3}' hybrid.all.asm.p_ctg.gfa > {output}
        """

rule hifiasm_single:
    input:
        hifi=W+"fastp/single.fq",
        ont=W+"fastp/ont.fq"
    output:
        W+"hifiasm_single/single.all.asm.p_ctg.fa"
    params:
        W + "hifiasm_single"
    shell:
        """
        cd {params}
        hifiasm -o single.all.asm --primary -t 96 --ul {input.ont} -k 63 {input.hifi}
        awk '/^S/{print ">"$2;print $3}' single.all.asm.p_ctg.gfa > {output}
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

rule rm_mt_cp_hifiasm_hybrid:
    input:
        hybrid=W+"hifiasm_hybrid/hybrid.all.asm.p_ctg.fa",
        mt=mitochondrion,
        cp=chloroplast
    output:
        W+"hifiasm_hybrid/hybrid.remove_cp_mt.fa"
    params:
        W+"hifiasm_hybrid"
    shell:
        """
        cd {params}
        minimap2 -t 96 -x asm5 {input.mt} {input.hybrid}> mitochondrion.paf
        minimap2 -t 96 -x asm5 {input.cp} {input.hybrid}> chloroplast.paf
        python gemma_los.py mitochondrion.paf > mitochondrion.txt
        python gemma_los.py chloroplast.paf > chloroplast.txt
        seqkit grep -v -f chloroplast.txt {input.hybrid} > wheat_remove_cp.fa
        seqkit grep -v -f mitochondrion.txt wheat_remove_cp.fa > {output}
        """

rule rm_mt_cp_hifiasm_single:
    input:
        single=W+"hifiasm_single/single.all.asm.p_ctg.fa",
        mt=mitochondrion,
        cp=chloroplast
    output:
        W+"hifiasm_single/single.remove_cp_mt.fa"
    params:
        W+"hifiasm_single"
    shell:
        """
        cd {params}
        minimap2 -t 96 -x asm5 {input.mt} {input.single}> mitochondrion.paf
        minimap2 -t 96 -x asm5 {input.cp} {input.single}> chloroplast.paf
        python gemma_los.py mitochondrion.paf > mitochondrion.txt
        python gemma_los.py chloroplast.paf > chloroplast.txt
        seqkit grep -v -f chloroplast.txt {input.single} > wheat_remove_cp.fa
        seqkit grep -v -f mitochondrion.txt wheat_remove_cp.fa > {output}
        """

rule hicpro_hybrid:
    input:
        hic=hic_hybrid_dir,
        ref=W+"hifiasm_hybrid/hybrid.remove_cp_mt.fa"
    output:
        W+"hic_hybrid/hic_hybrid.bam"
    params:
        dir=W+"hic_hybrid",
        prefix="hybrid.remove_cp_mt",
        spart_dir=SPART_dir,
        enzyme=hic_hybrid_enzyme
    shell:
        """
        cd {params.dir}
        ln -s {input.ref} ./
        bowtie2-build --large-index --threads 96 {params.prefix}.fa {params.prefix}
        samtools faidx {params.prefix}.fa
        awk '{print $1 "\t" $2}' {params.prefix}.fa.fai > genome_sizes.bed
        python ./HiC-Pro/bin/utils/digest_genome.py -r ^{params.enzyme} -o enzyme.bed {params.prefix}.fa
        makeblastdb -in {params.prefix}.fa -dbtype nucl -parse_seqids -out {params.prefix}
        cp {params.spart_dir}/01_Contig_scaffolding/hicpro_config.txt ./
        sed -i 's#^N_CPU = #N_CPU = 96#g' hicpro_config.txt
        sed -i 's#^BOWTIE2_IDX_PATH = #BOWTIE2_IDX_PATH = {params.dir}#g' hicpro_config.txt
        sed -i 's#^REFERENCE_GENOME = #REFERENCE_GENOME = {params.prefix}#g' hicpro_config.txt
        sed -i 's#^GENOME_SIZE = #GENOME_SIZE = {params.dir}/genome_sizes.bed#g' hicpro_config.txt
        sed -i 's#^GENOME_FRAGMENT = #GENOME_FRAGMENT = {params.dir}/enzyme.bed#g' hicpro_config.txt
        HiC-Pro -i {input.hic} -c hicpro_config.txt -o {params.dir}
        cd bowtie_results/bwt2
        for item in dir {params.dir}/bowtie_results/bwt2/*/*.bwt2pairs.bam; do samtools sort -m 1500M -n -@ 96 $item > $item.bam; done
        samtools merge -@ 96 -o {output} {params.dir}/bowtie_results/bwt2/*/*.bwt2pairs.bam.bam
        """

rule hicpro_single:
    input:
        hic=hic_single_dir,
        ref=W+"hifiasm_single/single.remove_cp_mt.fa"
    output:
        W+"hic_single/hic_single.bam"
    params:
        dir=W+"hic_single",
        prefix="single.remove_cp_mt",
        spart_dir=SPART_dir,
        enzyme=hic_single_enzyme
    shell:
        """
        cd {params.dir}
        ln -s {input.ref} ./
        bowtie2-build --large-index --threads 96 {params.prefix}.fa {params.prefix}
        samtools faidx {params.prefix}.fa
        awk '{print $1 "\t" $2}' {params.prefix}.fa.fai > genome_sizes.bed
        python ./HiC-Pro/bin/utils/digest_genome.py -r ^{params.enzyme} -o enzyme.bed {params.prefix}.fa
        makeblastdb -in {params.prefix}.fa -dbtype nucl -parse_seqids -out {params.prefix}
        cp {params.spart_dir}/01_Contig_scaffolding/hicpro_config.txt ./
        sed -i 's#^N_CPU = #N_CPU = 96#g' hicpro_config.txt
        sed -i 's#^BOWTIE2_IDX_PATH = #BOWTIE2_IDX_PATH = {params.dir}#g' hicpro_config.txt
        sed -i 's#^REFERENCE_GENOME = #REFERENCE_GENOME = {params.prefix}#g' hicpro_config.txt
        sed -i 's#^GENOME_SIZE = #GENOME_SIZE = {params.dir}/genome_sizes.bed#g' hicpro_config.txt
        sed -i 's#^GENOME_FRAGMENT = #GENOME_FRAGMENT = {params.dir}/enzyme.bed#g' hicpro_config.txt
        HiC-Pro -i {input.hic} -c hicpro_config.txt -o {params.dir}
        cd bowtie_results/bwt2
        for item in dir {params.dir}/bowtie_results/bwt2/*/*.bwt2pairs.bam; do samtools sort -m 1500M -n -@ 96 $item > $item.bam; done
        samtools merge -@ 96 -o {output} {params.dir}/bowtie_results/bwt2/*/*.bwt2pairs.bam.bam
        """

rule yahs_hybrid:
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
        yahs -e {params.enzyme} {input.ref} {input.bam} -o {params.prefix}
        cp {params.dir}/yahs_bam_scaffolds_final.fa ./yahs_hybrid.fa
        """

rule yahs_single:
    input:
        bam=W+"hic_single/hic_single.bam",
        ref=W+"hifiasm_single/single.remove_cp_mt.fa"
    output:
        W + "yahs_single/yahs_single.fa"
    params:
        dir = W + "yahs_single",
        prefix = "single_bam",
        enzyme = hic_single_enzyme
    shell:
        """
        cd {params.dir}
        yahs -e {params.enzyme} {input.ref} {input.bam} -o {params.prefix}
        cp {params.dir}/yahs_bam_scaffolds_final.fa ./yahs_single.fa
        """

rule ragtag_patch_single_hybrid:
    input:
        single=W + "yahs_single/yahs_single.fa",
        hybrid=W + "yahs_hybrid/yahs_hybrid.fa"
    output:
        W + "patch_hybrid/patch_single_hybrid.fa"
    params:
        dir = W + "patch_hybrid",
        prefix = "single_hybrid",
    shell:
        """
        cd {params.dir}
        wfmash {input.single} {input.hybrid} > {params.prefix}.paf 
        mkdir ragtag_output
        cd ragtag_output
        ln -s ../{params.prefix}.paf ragtag.patch.asm.paf
        cd ..
        ragtag.py patch -f 10000 --remove-small {input.single} {input.hybrid}
        cp {params.dir}/ragtag_output/ragtag.patch.fasta {output}
        """

rule ragtag_patch_flye:
    input:
        single_hybrid=W + "patch_hybrid/patch_single_hybrid.fa",
        flye=W + "flye/assembly.fasta"
    output:
        W + "patch_flye/patch_single_hybrid_flye.fa"
    params:
        dir = W + "patch_flye",
        prefix = "single_hybrid_flye",
    shell:
        """
        cd {params.dir}
        wfmash {input.single_hybrid} {input.flye} > {params.prefix}.paf 
        mkdir ragtag_output
        cd ragtag_output
        ln -s ../{params.prefix}.paf ragtag.patch.asm.paf
        cd ..
        ragtag.py patch -f 10000 --remove-small {input.single_hybrid} {input.flye}
        cp {params.dir}/ragtag_output/ragtag.patch.fasta {output}
        """

rule ragtag_patch_verkko
    input:
        single_hybrid_flye=W + "patch_hybrid/patch_single_hybrid.fa",
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
        wfmash {input.single_hybrid_flye} {input.verkko} > {params.prefix}.paf 
        mkdir ragtag_output
        cd ragtag_output
        ln -s ../{params.prefix}.paf ragtag.patch.asm.paf
        cd ..
        ragtag.py patch -f 10000 --remove-small {input.single_hybrid_flye} {input.verkko}
        cp {params.dir}/ragtag_output/ragtag.patch.fasta {output.ref}
        bwa-mem2.avx512bw index {input.ref}
        meryl count k=27 output merylDB {output.ref}
        meryl print greater-than distinct=0.9998 merylDB > {output.txt}
        """

rule minimap_cu_mix:
    input:
        fq=W+"hifi_mix_reads/{hifi_mix}_q40l15k.fastq",
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

rule minimap_cu_single:
    input:
        fq=W+"hifi_single_reads/{hifi_single}_q40l15k.fastq",
        ref = W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        txt = W + "repetitive_k27.txt"
    output:
        W+"hifi_single_winnowmap/{hifi_single}_q40l15k.sam"
    benchmark:
        W+"benchmarks/hifi_single_winnowmap/{hifi_single}.benchmark.txt"
    shell:
        # meryl count k=27 output merylDB CS_ISSA.fasta
        #meryl print greater-than distinct=0.9998 merylDB > repetitive_k27.txt
        "winnowmap --MD -W {input.txt} -ax map-pb -H -K 1500M -k 27 -w27 -t32 {input.ref} {input.fq} > {output}"

rule hifi_mix_sort:
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

rule filter:
    input:
        W+"hifi_mix_sort/{hifi_mix}_q40l15k.bam"
    output:
        W+"hifi_mix_sort_filter/{hifi_mix}_q40l15k.bam"
    benchmark:
        W + "benchmarks/hifi_mix_sort_filter/{hifi_mix}.benchmark.txt"
    shell:
        "samtools view -@32 -F0x104 -hb {input} > {output}"

rule filter_merge:
    input:
        expand(W+"hifi_mix_sort_filter/{hifi_mix}_q40l15k.bam",hifi_mix=hifi_mix)
    output:
        W+"hybrid/hybrid.bam"
    benchmark:
        W + "benchmarks/hybrid/hybrid.benchmark.txt"
    shell:
        "samtools merge -@ 128 -l 0 {output} {input}"

rule pcr_free_hybrid:
    input:
        fa=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        r1=pcrfree_hybrid_r1,
        r2=pcrfree_hybrid_r2
    output:
        W+"hybrid_hifi_pcr/pcr.bam"
    shell:
       "bwa-mem2.avx512bw mem -t 96 {input.fa} {input.r1} {input.r2}|samtools view -@ 96 -b -|samtools sort -@ 96 -m 30G -o {output} -"

rule pcr_free_single:
    input:
        fa=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        r1=pcrfree_single_r1,
        r2=pcrfree_single_r2
    output:
        W+"single_hifi_pcr/pcr.bam"
    shell:
       "bwa-mem2.avx512bw mem -t 96 {input.fa} {input.r1} {input.r2}|samtools view -@ 96 -b -|samtools sort -@ 96 -m 30G -o {output} -"

rule filter_merge_hybrid:
    input:
        hifi=expand(W+"hifi_mix_sort_filter/{hifi_mix}_q40l15k.bam",hifi_mix=hifi_mix),
        pcr=W+"hybrid_hifi_pcr/pcr.bam"
    output:
        W+"hybrid_hifi_pcr/hybrid.bam"
    benchmark:
        W + "benchmarks/hybrid_pcr/hybrid.benchmark.txt"
    shell:
        "samtools merge -@ 128 -l 0 {output} {input.hifi} {input.pcr}"

rule hifi_single_sort:
    input:
        W+"hifi_single_winnowmap/{hifi_single}_q40l15k.sam"
    output:
        W+"hifi_single_sort/{hifi_single}_q40l15k.bam"
    params:
        W + "patch_verkko/patch_single_hybrid_flye_verkko.fa.fai"
    benchmark:
        W + "benchmarks/hifi_single_sort/{hifi_single}.benchmark.txt"
    shell:
        "samtools view -@32 -bt {params} {input}|samtools sort -@32 -m1500M -O bam -o {output} -"

rule filter_single:
    input:
        W+"hifi_single_sort/{hifi_single}_q40l15k.bam"
    output:
        W+"hifi_single_sort_filter/{hifi_single}_q40l15k.bam"
    benchmark:
        W + "benchmarks/hifi_single_sort_filter/{hifi_single}.benchmark.txt"
    shell:
        "samtools view -@ 32 -F0x104 -hb {input} > {output}"

rule filter_merge_single:
    input:
        expand(W+"hifi_single_sort_filter/{hifi_single}_q40l15k.bam",hifi_single=d)
    output:
        W+"single/single.bam"
    benchmark:
        W + "benchmarks/single/single.benchmark.txt"
    shell:
        "samtools merge -@ 128 -l 0 {output} {input}"

rule filter_merge_single_pcr:
    input:
        hifi=expand(W+"hifi_single_sort_filter/{hifi_single}_q40l15k.bam",hifi_single=d),
        pcr=W+"single_hifi_pcr/pcr.bam"
    output:
        W+"single_hifi_pcr/hybrid.bam"
    benchmark:
        W + "benchmarks/single_pcr/hybrid.benchmark.txt"
    shell:
        "samtools merge -@ 128 -l 0 {output} {input.hifi} {input.pcr}"


rule minimap_cu_ont:
    input:
        fq=W+"ont/reads/{e}_q10l120k.fastq",
        ref=W + "patch_verkko/patch_single_hybrid_flye_verkko.fa",
        txt=W+"repetitive_k27.txt"
    output:
        W+"ont_winnowmap/{e}/{e}_q10l120k.sam"
    benchmark:
        W+"benchmarks/ont_winnowmap/{e}.benchmark.txt"
    shell:
        "winnowmap --MD -W {input.txt} -ax map-ont -H -K 1500M -k 27 -w27 -t32 {input.ref} {input.fq} > {output}"

rule sort:
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

rule filter_ont:
    input:
        W+"ont_sort/{e}/{e}_q10l120k.bam"
    output:
        W+"ont_filter/{e}_q10l120k.bam"
    benchmark:
        W + "benchmarks/ont_filter/{e}.benchmark.txt"
    shell:
        "samtools view -@ 128 -F0x104 -hb {input} > {output}"

rule merge:
    input:
        expand(W+"ont_filter/{e}_q10l120k.bam",e=e)
    output:
        W + "ont_merge/q10l120k.bam"
    benchmark:
        W + "benchmarks/ont_merge/benchmark.txt"
    shell:
        "samtools merge -@ 128 -o {output} {input}"