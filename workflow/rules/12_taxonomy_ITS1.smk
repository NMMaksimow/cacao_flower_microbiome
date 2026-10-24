# =============================================================================
# 12. ITS1 taxonomy: ASVs and OTUs, each with both UNITE classifiers
# =============================================================================
# As v1 (qiime2/scripts/08b for ASVs, rule taxonomy_its1_otu for OTUs): default
# confidence 0.7 and one worker (--p-n-jobs 1): each worker loads its own copy of
# the classifier. Unlike v1, --p-read-orientation same instead of the default
# "auto", which guesses the orientation from a subset of reads (see rule 09); the
# ITS1 reads start right after the ITS1F primer, the orientation of the UNITE
# references. Both classifiers are compared in QC 07; the main one is set in the
# config (taxonomy_main) and exported in rule 14.

ITS1_FEATURES = {"asv": "CFM_ITS1_dada2_repseqs", "otu97": "CFM_ITS1_otu97_repseqs"}


rule classify_ITS1:
    input:
        classifier=lambda wildcards: UNITE[wildcards.classifier],
        reads=lambda wildcards: f"{OUT_DIR}/denoise/{ITS1_FEATURES[wildcards.features]}.qza",
    output:
        qza=f"{OUT_DIR}/taxonomy/CFM_ITS1_{{features}}_taxonomy_{{classifier}}.qza",
        tsv=f"{OUT_DIR}/taxonomy/CFM_ITS1_{{features}}_taxonomy_{{classifier}}.tsv",
    wildcard_constraints:
        features="asv|otu97",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: 1
    resources:
        partition="bigmem",   # s_all on the 12,388 ASVs ran out of memory at 150 GB
        mem_mb=300000,
        runtime=480,
        tmp_gb=20,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.reads} \
            --p-read-orientation same \
            --p-n-jobs 1 \
            --o-classification {output.qza}
        export_dir=$(mktemp -d)
        {params.r_isolation} qiime tools export --input-path {output.qza} --output-path $export_dir
        mv $export_dir/taxonomy.tsv {output.tsv}
        rm -r $export_dir
        """
