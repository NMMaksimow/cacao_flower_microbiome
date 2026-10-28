# =============================================================================
# QC 06. Comparison of the four SILVA 144 classifiers (choice of 16S classifier)
# =============================================================================
# Side check; nothing downstream uses its output.


rule qc_silva144_comparison:
    input:
        table=f"{OUT_DIR}/05_merged/{config['main_rung']}/16S/feature-table.biom",
        asvs=f"{OUT_DIR}/05_merged/{config['main_rung']}/16S/dna-sequences.fasta",
        expected=f"{workflow.basedir}/resources/mock/{EXPECTED_MOCK['16S']}",
        taxonomies=expand(f"{OUT_DIR}/taxonomy/CFM_16S_taxonomy_{{classifier}}.tsv",
                          classifier=list(SILVA144)),
    output:
        summary=f"{OUT_DIR}/qc/06_silva144_comparison.tsv",
        details=f"{OUT_DIR}/qc/06_silva144_details.tsv",
        figure=f"{OUT_DIR}/qc/06_silva144_comparison.png",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        script=f"{workflow.basedir}/scripts/qc_classifier_comparison.py",
        min_mock=config["dada2"]["mock_min_pairs"],
        named=" ".join(f"{c}={OUT_DIR}/taxonomy/CFM_16S_taxonomy_{c}.tsv" for c in SILVA144),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        python {params.script} 16S {input.table} {input.asvs} {input.expected} {params.min_mock} \
            {output.summary} {output.details} {output.figure} {params.named}
        """
