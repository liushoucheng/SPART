# Install
git clone https://github.com/liushoucheng/SPART.git

cd SPART

conda env create -f SPART.yaml

conda activate spart

# Dependencies

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

# Running pipeline with snakemake
Exclude Verkko,Bionano DLS Map,Telomere determination and patch,Centromeric region analysis,Variant calls and Evaluation

sed -i "s#^ SPART_PATH# ${PWD}#g" conf_ck.yaml

HiC_enzyme=" GATC"

sed -i "s#^ hic_sca_enzyme# ${HiC_enzyme}#g" conf_ck.yaml

HiC_ligation_site=" GATCGATC"

sed -i "s#^ hic_sca_ligation_site# ${HiC_ligation_site}#g" conf_ck.yaml

snakemake -s SPART.py --cluster-config clust.json --configfile conf_ck.yaml --cluster '{cluster.account}' --jobs $threads --rerun-incomplete --restart-times 1 -np --rulegraph |dot -Tpng > rule.png

configfile:The config file can be used to define a dictionary of configuration parameters and their values.

cluster-config:A JSON or YAML file that defines the wildcards used in 'cluster'for specific rules.
<div align=center>
<img src="https://github.com/liushoucheng/SPART/blob/main/pic/rule.png">
</div>

## Output files
please see the complete [documentation]( https://github.com/liushoucheng/SPART/tree/main/exmple).

