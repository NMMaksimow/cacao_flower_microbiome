# =============================================================================
# QC 07. Comparison of the two UNITE classifiers on the ITS1 OTUs (choice of ITS1
# classifier; decision rule fixed in docs/v2_rerun_plan.md before the results)
# =============================================================================
# Side check; nothing downstream uses its output. Same script and metrics as QC 06;
# the reference classifier is s_all (v1's).


rule qc_unite_comparison:
    input:
        table=f"{EXPORT_DIR}/CFM_ITS1_otu97_table/feature-table.biom",
        otus=f"{EXPORT_DIR}/CFM_ITS1_otu97_repseqs/dna-sequences.fasta",
        expected=f"{workflow.basedir}/resources/mock/{EXPECTED_MOCK['ITS1']}",
        taxonomies=expand(f"{OUT_DIR}/taxonomy/CFM_ITS1_otu97_taxonomy_{{classifier}}.tsv",
                          classifier=list(UNITE)),
    output:
        summary=f"{OUT_DIR}/qc/07_unite_comparison.tsv",
        details=f"{OUT_DIR}/qc/07_unite_details.tsv",
        figure=f"{OUT_DIR}/qc/07_unite_comparison.png",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        script=f"{workflow.basedir}/scripts/qc_classifier_comparison.py",
        min_mock=config["dada2"]["mock_min_pairs"],
        named=" ".join(f"{c}={OUT_DIR}/taxonomy/CFM_ITS1_otu97_taxonomy_{c}.tsv" for c in UNITE),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        python {params.script} ITS1 {input.table} {input.otus} {input.expected} {params.min_mock} \
            {output.summary} {output.details} {output.figure} {params.named}
        """
