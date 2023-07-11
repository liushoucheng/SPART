import os
import re

b={}
d={}
e={}
W=config["WORKDIR"]

rule final:
    input:
        W + "output/SV_SNV/merfin_sv_snv_consensus.fasta"

rule sv_sniffles_hybrid:
    input:
        W+"hybrid/hybrid.bam"
    output:
        W + "output/SV/sniffles_hybrid.vcf"
    shell:
        "sniffles --threads 128 --input {input} --vcf {output}"

rule sv_sniffles_single:
    input:
        W+"single/single.bam"
    output:
        W + "output/SV/sniffles_single.vcf"
    shell:
        "sniffles --threads 128 --input {input} --vcf {output}"

rule sv_sniffles_ont:
    input:
        W + "ont_merge/q10l120k.bam"
    output:
        W + "output/SV/sniffles_ont.vcf"
    shell:
        "sniffles --threads 128 --input {input} --vcf {output}"

#################call SV
rule sv_cutesv_hybrid:
    input:
        bam=W+"hybrid/hybrid.bam",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SV/cutesv_hybrid.vcf"
    params:
        W + "output/SV/cutesv_hybrid"
    shell:
        "cuteSV --threads 128 {input.bam} {input.ref} {output} {params} --genotype"

rule sv_cutesv_hybrid_suggest:
    input:
        bam=W+"hybrid/hybrid.bam",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SV/cutesv_hybrid_suggest.vcf"
    params:
        W + "output/SV/cutesv_hybrid_suggest"
    shell:
        "cuteSV --threads 128 {input.bam} {input.ref} {output} {params} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype"

rule sv_cutesv_single:
    input:
        bam=W+"single/single.bam",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SV/cutesv_single.vcf"
    params:
        W + "output/SV/cutesv_single"
    shell:
        "cuteSV --threads 128 {input.bam} {input.ref} {output} {params} --genotype"

rule sv_cutesv_single_suggest:
    input:
        bam=W+"single/single.bam",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SV/cutesv_single_suggest.vcf"
    params:
        W + "output/SV/cutesv_single_suggest"
    shell:
        "cuteSV --threads 128 {input.bam} {input.ref} {output} {params} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype"

rule sv_cutesv_ont:
    input:
        bam=W + "ont_merge/q10l120k.bam",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SV/cutesv_ont.vcf"
    params:
        W + "output/SV/cutesv_ont"
    shell:
        "cuteSV --threads 128 {input.bam} {input.ref} {output} {params} --genotype"

rule sv_cutesv_ont_suggest:
    input:
        bam=W + "ont_merge/q10l120k.bam",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SV/cutesv_ont_suggest.vcf"
    params:
        W + "output/SV/cutesv_ont_suggest"
    shell:
        "cuteSV --threads 128 {input.bam} {input.ref} {output} {params} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --genotype"

rule hifi_ls:
    input:
        hifi1=W + "output/SV/sniffles_hybrid.vcf",
        hifi2=W + "output/SV/cutesv_hybrid.vcf",
        hifi3=W + "output/SV/cutesv_hybrid_suggest.vcf",
        hifi4=W + "output/SV/sniffles_single.vcf",
        hifi5=W + "output/SV/cutesv_single.vcf",
        hifi6=W + "output/SV/cutesv_single_suggest.vcf"
    output:
        W + "output/SV/hifi_ls.txt"
    shell:
        "ls {input.hifi1} {input.hifi2} > {output}"

rule ont_ls:
    input:
        ont1=W + "output/SV/{chr}/{chr}_q10l120k_cutesv_ont.vcf",
        ont2=W + "output/SV/{chr}/{chr}_q10l120k_cutesv_ont_suggest.vcf",
        ont3=W + "output/SV/sniffles_ont.vcf"
    output:
        W + "output/SV/ont_ls.txt"
    shell:
        "ls {input.ont1} {input.ont2} > {output}"

rule jasmine_hifi:
    input:
        W + "output/SV/hifi_ls.txt"
    output:
        W + "output/SV/jasmine_hifi.vcf"
    shell:
        "jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 threads=128 min_support=1 --output_genotypes file_list={input} out_file={output}"

rule jasmine_ont:
    input:
        W + "output/SV/ont_ls.txt"
    output:
        W + "output/SV/jasmine_ont.vcf"
    shell:
        "jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 threads=128 min_support=1 --output_genotypes file_list={input} out_file={output}"

rule ls_hifi_ont:
    input:
        hifi=W + "output/SV/jasmine_hifi.vcf",
        ont=W + "output/SV/jasmine_ont.vcf"
    output:
        W + "output/hifi_ont_ls.txt"
    shell:
        "ls {input.hifi} {input.ont} > {output}"

rule jasmine_hifi_ont:
    input:
        W + "output/SV/hifi_ont_ls.txt"
    output:
        W + "output/SV/jasmine_hifi_ont.vcf"
    shell:
        "jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 threads=96 min_support=2 --output_genotypes file_list={input} out_file={output}"

############################call SNV

rule vcf_filter_ont:
    input:
        W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz"
    output:
        W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.PASS.gq25.gt10.vcf.gz"
    shell:
        "bcftools view -f ""PASS"" -e 'FORMAT/GQ<=25 | FORMAT/DP<=10' -Oz {input} > {output}"

rule vcf_filter_ont_index:
    input:
        W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.PASS.gq25.gt10.vcf.gz"
    output:
        W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.PASS.gq25.gt10.vcf.gz.csi"
    shell:
        "bcftools index -c {input}"

rule vcf_filter_hifi:
    input:
        W + "output/hybrid.vcf"
    output:
        W + "output/hybrid.PASS.gq25.gt10.vcf.gz"
    shell:
        "bcftools view -f ""PASS"" -e 'FORMAT/GQ<=25 | FORMAT/DP<=10' -Oz {input} > {output}"

rule vcf_filter_hifi_index:
    input:
        W + "output/hybrid.PASS.gq25.gt10.vcf.gz"
    output:
        W + "output/hybrid.PASS.gq25.gt10.vcf.gz.csi"
    shell:
        "bcftools index -c {input}"

rule hap:
    input:
        hifi=W + "output/hybrid.PASS.gq25.gt10.vcf.gz",
        hifi_csi=W + "output/hybrid.PASS.gq25.gt10.vcf.gz.csi",
        ont=W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.PASS.gq25.gt10.vcf.gz",
        ont_csi=W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.PASS.gq25.gt10.vcf.gz.csi",
        ref="CS_ISSA.fasta"
    output:
        W + "output/SNV/HAPPY.vcf.gz",
    params:
         W + "output/SNV/HAPPY"
    shell:
        "python hap.py {input.hifi} {input.ont} -r {input.ref} -o {params} --pass-only --threads 128"

rule vcf_merge_t2t:
    input:
        hifi=W + "output/hybrid.PASS.gq25.gt10.vcf.gz",
        ont=W + "output/pepper_deepvariant_output/intermediate_files/PEPPER_VARIANT_FULL.PASS.gq25.gt10.vcf.gz",
        hap=W + "output/SNV/HAPPY.vcf.gz",
    output:
        W + "output/SNV/MERGED_SMALL_VARIANTS.vcf.gz"
    shell:
        "python3 vcf_merge_t2t.py -v1 {input.hifi} -v2 {input.ont} -hv {input.hap} -o {output}"

rule gunzip:
    input:
        W + "output/SNV/MERGED_SMALL_VARIANTS.vcf.gz"
    output:
        W + "output/SNV/MERGED_SMALL_VARIANTS.vcf"
    shell:
        "gunzip -d -c {input} > {output}"

rule meryl_count:
    input:
        "CS_ISSA.fasta"
    output:
        directory(W+"meryl/merylDB_k21"),
    params:
        k="21",
        dir=W+"meryl/merylDB{chr}_k21"
    threads: 128
    shell:
        "/home/liusc/software/meryl-1.4/bin/meryl count k={params.k} threads={threads} {input} output {params.dir}"

rule merfin_snv:
    input:
        ref="CS_ISSA.fasta",
        seqmers=W+"meryl/merylDB_k21",
        vcf=W + "output/SNV/MERGED_SMALL_VARIANTS.vcf"
    output:
        W + "output/SNV/merfin_snv.filter.vcf"
    params:
        W + "output/SNV/merfin_snv"
    shell:
        "merfin -strict -threads 128 -sequence {input.ref} -seqmers {input.seqmers} -readmers single.hifi40_cspcrfree.k21.gt1.meryl -peak 88.3 -prob lookup_table.txt -vcf {input.vcf} -output {params}"

rule merfin_sv:
    input:
        ref="CS_ISSA.fasta",
        seqmers=W+"meryl/merylDB_k21",
        vcf=W + "output/SV/jasmine_hifi_ont.vcf"
    output:
        W + "output/SV/merfin_sv.filter.vcf"
    params:
        W + "output/SV/merfin_sv"
    shell:
        "merfin -strict -threads 128 -sequence {input.ref} -seqmers {input.seqmers} -readmers single.hifi40_cspcrfree.k21.gt1.meryl -peak 88.3 -prob single.hifi40_cspcrfree.k21/lookup_table.txt -vcf {input.vcf} -output {params}"

rule cut:
    input:
        W + "output/SV/merfin_sv.filter.vcf"
    output:
        W + "output/SV/merfin_sv10.filter.vcf"
    shell:
        "cut -f 1-10 {input} > {output}"

rule ls_sv_snv:
    input:
        snv=W + "output/SNV/merfin_snv.filter.vcf",
        sv=W + "output/SV/merfin_sv10.filter.vcf"
    output:
        W + "output/SV_SNV/lst.txt"
    shell:
        "ls {input.snv} {input.sv} > {output}"
rule jasmine_merge:
    input:
        W + "output/SV_SNV/lst.txt"
    output:
        W + "output/SV_SNV/SV_SNV.vcf"
    shell:
        "jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 threads=96 --output_genotypes file_list={input} out_file={output}"
rule merfin_sv_snv:
    input:
        ref="CS_ISSA.fasta",
        seqmers=W+"meryl/merylDB_k21",
        vcf=W + "output/SV_SNV/SV_SNV.vcf"
    output:
        W + "output/SV_SNV/merfin_sv_snv.filter.vcf"
    params:
        W + "output/SV_SNV/merfin_sv_snv"
    shell:
        "merfin -strict -threads 128 -sequence {input.ref} -seqmers {input.seqmers} -readmers single.hifi40_cspcrfree.k21.gt1.meryl -peak 88.3 -prob single.hifi40_cspcrfree.k21/lookup_table.txt -vcf {input.vcf} -output {params}"

rule cut_merfin_sv_snv:
    input:
        W + "output/SV_SNV/merfin_sv_snv.filter.vcf"
    output:
        W + "output/SV_SNV/merfin_sv_snv.filter10.vcf"
    shell:
        "cut -f 1-10 {input} > {output}"

rule view_merfin_sv_snv:
    input:
        W + "output/SV_SNV/merfin_sv_snv.filter10.vcf"
    output:
        W + "output/SV_SNV/merfin_sv_snv.filter.vcf.gz"
    shell:
        "bcftools view -Oz {input} > {output}"

rule view_merfin_sv_snv_index:
    input:
        W + "output/SV_SNV/merfin_sv_snv.filter.vcf.gz"
    output:
        W + "output/SV_SNV/merfin_sv_snv.filter.vcf.gz.csi"
    shell:
        "bcftools index -c {input}"

rule fasta:
    input:
        vcf=W + "output/SV_SNV/merfin_sv_snv.filter.vcf.gz",
        index=W + "output/SV_SNV/merfin_sv_snv.filter.vcf.gz.csi",
        ref="CS_ISSA.fasta",
    output:
        W + "output/SV_SNV/merfin_sv_snv_consensus.fasta"
    shell:
        "bcftools consensus -f {input.ref} -H 1 {input.vcf} > {output}"


#snakemake -s callsv_snv.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs 128 --rerun-incomplete --restart-times 1 -np