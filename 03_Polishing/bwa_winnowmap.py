import os
import re

b={}
hifi_single={}
hifi_mix={}
e={}
d={}
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



rule minimap_cu_mix:
    input:
        fq=W+"hifi_mix_reads/{hifi_mix}_q40l15k.fastq",
        ref="CS_ISSA.fasta",
        txt="repetitive_k27.txt"
    output:
        W+"hifi_mix_winnowmap/{hifi_mix}_q40l15k.sam"
    benchmark:
        W+"benchmarks/hifi_mix_winnowmap/{hifi_mix}.benchmark.txt"
    shell:
        # meryl count k=27 output merylDB CS_ISSA.fasta
        #meryl print greater-than distinct=0.9998 merylDB > repetitive_k27.txt
        "winnowmap --MD -W {input.txt} -ax map-pb -H -K 1500M -k 27 -w27 -t32 {input.ref} {input.fq} > {output}"

rule minimap_cu_single:
    input:
        fq=W+"hifi_single_reads/{hifi_single}_q40l15k.fastq",
        ref = "CS_ISSA.fasta",
        txt = "repetitive_k27.txt"
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
        "CS_ISSA.fasta.fai"
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


#"/home/liusc/lxp/software/bwa-mem2/bwa-mem2.avx512bw index {input}"

rule pcr_free_hybrid:
    input:
        fa="CS_ISSA.fasta",
        r1="hybrid_1.clean.fq.gz",
        r2="hybrid_2.clean.fq.gz"
    output:
        W+"hybrid_hifi_pcr/pcr.bam"
    shell:
       "bwa-mem2.avx512bw mem -t 96 {input.fa} {input.r1} {input.r2}|samtools view -@ 96 -b -|samtools sort -@ 96 -m 30G -o {output} -"

rule pcr_free_single:
    input:
        fa="CS_ISSA.fasta",
        r1="single_1.clean.fq.gz",
        r2="single_2.clean.fq.gz"
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
        "/data/liusc/lixp/wheat/result/ref/m219/m219.fasta.fai"
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
        ref="CS_ISSA.fasta",
        txt="repetitive_k27.txt"
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
        "CS_ISSA.fasta.fai"
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

# snakemake -s bwa_winnowmap.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs 40 --rerun-incomplete --restart-times 1 -np