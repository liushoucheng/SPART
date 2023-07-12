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
W=config["WORKDIR"]
DIR=config["DIR"]
DIRs=config["DIRs"]
DIRont=config["DIRont"]

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
rule hybrid_fastp:
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

    params:
        dir=W+"hic_hybrid"
        prefix="hybrid.remove_cp_mt"
        spart_dir=SPART_dir
    shell:
        """
        cd {params.dir}
        ln -s {input.ref} ./
        bowtie2-build --large-index --threads 96 {params.prefix}.fa {params.prefix}
        samtools faidx {params.prefix}.fa
        awk '{print $1 "\t" $2}' {params.prefix}.fa.fai > genome_sizes.bed
        python ./HiC-Pro/bin/utils/digest_genome.py -r ^GATC -o wheat_DpnII.bed {params.prefix}.fa
        makeblastdb -in {params.prefix}.fa -dbtype nucl -parse_seqids -out {params.prefix}
        cp {params.spart_dir}/01_Contig_scaffolding/hicpro_config.txt ./
        sed -i 's#^N_CPU = #N_CPU = 96#g' hicpro_config.txt
        sed -i 's#^protein=#protein={params.pro}#g' maker_opts.ctl
        sed -i 's#^rmlib=#rmlib={params.te}#g' maker_opts.ctl
        HiC-Pro -i {input.hic} -c $hicpro_config -o $hicpro_outdir
        """
