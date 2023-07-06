#!/bin/sh

threads=$1


python pipelineCL.py -Tn $threads -i 5 -b /home/liusc/proj/wheat/rawdata/bionano/output.bnx -r wheat.all.asm.p_ctg_50k_CTTAAG_0kb_0labels.cmap -l all_p50 -e w -a /home/liusc/proj/wheat/result/bionano/optArguments_nonhaplotype_noES_noCut_BG_DLE1_saphyr.xml -t /home/bionano/tools/pipeline/1.0/RefAligner/1.0 -y -z --species-reference other -C /home/bionano/tools/pipeline/1.0/Pipeline/1.0/clusterArguments_saphyr.xml -F 1
perl /home/bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold.pl -n /home/liusc/proj/wheat/result/hifiasm/hifi.20.30x/wheat.asm.p_ctg.fa -b /home/liusc/proj/wheat/result/bionano/all_p50/contigs/w_refineFinal1/W_REFINEFINAL1.cmap -u CTTAAG -c /home/bionano/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -r /home/bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -o bio_hybrid -B 2 -N 2 -f