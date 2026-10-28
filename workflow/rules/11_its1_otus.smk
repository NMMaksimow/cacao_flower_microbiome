# =============================================================================
# 11. ITS1 OTUs: vsearch de-novo clustering of the ASVs at 97% identity
# =============================================================================
# As v1 (workflow/Snakefile, rule cluster_its1_otus): one fungal strain can give
# several ITS1 ASVs (intragenomic length polymorphism), 97% OTUs collapse them.
# Other options default (strand plus). The OTU ID is the ID of its centroid ASV.
# Same artifact names and folder as v1 (qiime2/denoise/), under qiime2_v2/.


rule cluster_its1_otus:
    input:
        table=f"{OUT_DIR}/denoise/CFM_ITS1_dada2_table.qza",
        repseqs=f"{OUT_DIR}/denoise/CFM_ITS1_dada2_repseqs.qza",
    output:
        table=f"{OUT_DIR}/denoise/CFM_ITS1_otu97_table.qza",
        repseqs=f"{OUT_DIR}/denoise/CFM_ITS1_otu97_repseqs.qza",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
        tmp_gb=20,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime vsearch cluster-features-de-novo \
            --i-table {input.table} \
            --i-sequences {input.repseqs} \
            --p-perc-identity 0.97 \
            --p-threads {threads} \
            --o-clustered-table {output.table} \
            --o-clustered-sequences {output.repseqs}
        """
