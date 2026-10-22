# =============================================================================
# 08. QIIME 2 import of the final ASV table and sequences (main rung only)
# =============================================================================
# Same artifact names and folder as v1 (qiime2/denoise/), under qiime2_v2/.


rule qiime2_import:
    input:
        f"{OUT_DIR}/05_merged/{config['main_rung']}/{{marker}}/seqtab_nochim.tsv",
    output:
        biom=f"{OUT_DIR}/05_merged/{config['main_rung']}/{{marker}}/feature-table.biom",
        fasta=f"{OUT_DIR}/05_merged/{config['main_rung']}/{{marker}}/dna-sequences.fasta",
        table=f"{OUT_DIR}/denoise/CFM_{{marker}}_dada2_table.qza",
        repseqs=f"{OUT_DIR}/denoise/CFM_{{marker}}_dada2_repseqs.qza",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
        script=f"{workflow.basedir}/scripts/seqtab_to_biom.py",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60,
        tmp_gb=5,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        python {params.script} {input} {output.biom} {output.fasta}
        {params.r_isolation} qiime tools import \
            --type 'FeatureTable[Frequency]' \
            --input-format BIOMV210Format \
            --input-path {output.biom} \
            --output-path {output.table}
        {params.r_isolation} qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {output.fasta} \
            --output-path {output.repseqs}
        """
