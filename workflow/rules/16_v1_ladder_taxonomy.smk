# =============================================================================
# 16. Ablation ladder: v1 sequences classified with the v2 main classifiers
# =============================================================================
# Separates the classifier effect from the processing effect in QC 10 (top taxa):
#   v1 sequences + v1 classifier  ->  v1 sequences + v2 classifier  ->  v2 (rung C)
# Same commands as rules 09 (16S) and 12 (ITS1). The v1 artifacts in qiime2/denoise/
# are only read. ITS1 uses the v1 97% OTUs, the level of the main ITS1 data.


rule classify_v1_16S:
    input:
        classifier=SILVA144[MAIN_16S],
        reads="qiime2/denoise/CFM_16S_dada2_repseqs.qza",
    output:
        qza=f"{OUT_DIR}/taxonomy/CFM_16S_v1_taxonomy_{MAIN_16S}.qza",
        tsv=f"{OUT_DIR}/taxonomy/CFM_16S_v1_taxonomy_{MAIN_16S}.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
        n_jobs=config["classify_n_jobs"],
    threads: config["classify_n_jobs"]
    resources:
        partition="bigmem",   # as classify_16S
        mem_mb=300000,
        runtime=360,
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
            --p-n-jobs {params.n_jobs} \
            --o-classification {output.qza}
        export_dir=$(mktemp -d)
        {params.r_isolation} qiime tools export --input-path {output.qza} --output-path $export_dir
        mv $export_dir/taxonomy.tsv {output.tsv}
        rm -r $export_dir
        """


rule classify_v1_ITS1_otus:
    input:
        classifier=UNITE[MAIN_ITS1],
        reads="qiime2/denoise/CFM_ITS1_otu97_repseqs.qza",
    output:
        qza=f"{OUT_DIR}/taxonomy/CFM_ITS1_v1_otu97_taxonomy_{MAIN_ITS1}.qza",
        tsv=f"{OUT_DIR}/taxonomy/CFM_ITS1_v1_otu97_taxonomy_{MAIN_ITS1}.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: 1
    resources:
        partition="bigmem",   # as classify_ITS1
        mem_mb=150000,
        runtime=240,
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
