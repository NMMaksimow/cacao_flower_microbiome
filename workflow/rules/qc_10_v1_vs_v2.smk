# =============================================================================
# QC 10. v1 vs v2 comparison table with the ablation rungs (Phase 2.8)
# =============================================================================
# Read retention per step, ASV and OTU counts, primer remnants, mock composition,
# top genera; v1, rungs A/B/C and v1 sequences with the v2 classifiers (rule file
# 16). Side check; nothing downstream uses its output. v1 files are only read.


rule qc_v1_vs_v2:
    input:
        v1_stats=expand("qiime2/denoise/CFM_{marker}_dada2_stats.qza", marker=MARKERS),
        v1_tables=expand("qiime2/export/CFM_{marker}_dada2_table/feature-table.biom", marker=MARKERS),
        v1_its1_otus="qiime2/export/CFM_ITS1_otu97_table/feature-table.biom",
        steps=f"{OUT_DIR}/qc/03_primer_trimming_per_sample.tsv",
        remnants=f"{OUT_DIR}/qc/05_primer_remnants.tsv",
        seqtabs=expand(f"{OUT_DIR}/05_merged/{{rung}}/{{marker}}/seqtab_nochim.tsv", rung=list(RUNGS), marker=MARKERS),
        its1_otus=f"{EXPORT_DIR}/CFM_ITS1_otu97_table/feature-table.biom",
        taxonomy=[f"{MAIN_TAXONOMY_16S}.tsv",
                  f"{OUT_DIR}/taxonomy/CFM_ITS1_otu97_taxonomy_{MAIN_ITS1}.tsv",
                  f"{OUT_DIR}/taxonomy/CFM_16S_v1_taxonomy_{MAIN_16S}.tsv",
                  f"{OUT_DIR}/taxonomy/CFM_ITS1_v1_otu97_taxonomy_{MAIN_ITS1}.tsv"],
    output:
        table_16S=f"{OUT_DIR}/qc/10_v1_vs_v2_16S.tsv",
        table_ITS1=f"{OUT_DIR}/qc/10_v1_vs_v2_ITS1.tsv",
        figure=f"{OUT_DIR}/qc/10_v1_vs_v2.png",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        script=f"{workflow.basedir}/scripts/qc_v1_vs_v2.py",
        mock_dir=f"{workflow.basedir}/resources/mock",
        min_mock=config["dada2"]["mock_min_pairs"],
        rungs=",".join(f"{name}={'+'.join(r['lanes'])}:{r['errfun']}" for name, r in RUNGS.items()),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        python {params.script} {OUT_DIR} qiime2 {params.mock_dir} {params.min_mock} {params.rungs} \
            {MAIN_16S} {MAIN_ITS1} {output.table_16S} {output.table_ITS1} {output.figure}
        """
