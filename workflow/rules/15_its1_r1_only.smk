# =============================================================================
# 15. ITS1 sensitivity branch: forward reads (R1) only
# =============================================================================
# ITS1 amplicons longer than ~430 nt (e.g. Saccharomyces, ~440 nt) cannot merge with
# 2x250 reads, so the main ITS1 data (merged pairs) lose them. This branch denoises
# R1 alone for all samples, from the main pipeline's filtered reads and binned R1
# error models (nothing relearned), merges the runs, removes chimeras, clusters 97%
# OTUs and classifies with the main ITS1 classifier. QC 09 compares it with the
# merged-read ASVs; what is reported was fixed in advance (docs/v2_rerun_plan.md).
# The main ITS1 data stay the merged reads.

R1_DIR = f"{OUT_DIR}/r1_only"


rule dada2_denoise_r1:
    input:
        filtered=f"{OUT_DIR}/04_dada2/{{lane}}/ITS1/filtered",
        err_r1=f"{OUT_DIR}/04_dada2/{{lane}}/ITS1/errors_binned_R1.rds",
    output:
        seqtab=f"{R1_DIR}/{{lane}}/seqtab.rds",
        track=f"{R1_DIR}/{{lane}}/track.tsv",
    wildcard_constraints:
        lane="L1|L2",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        script=f"{workflow.basedir}/scripts/dada2_denoise_r1.R",
        max_length=config["r1_only"]["max_length"],
    threads: 16
    resources:
        mem_mb=32000,
        runtime=240,
        tmp_gb=5,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        {params.r_isolation} Rscript {params.script} \
            {input.filtered} {input.err_r1} {threads} {params.max_length} \
            {output.seqtab} {output.track}
        """


# Same script and settings as rule dada2_merge_chimeras (sum across runs, consensus)
rule dada2_merge_chimeras_r1:
    input:
        expand(f"{R1_DIR}/{{lane}}/seqtab.rds", lane=LANES),
    output:
        rds=f"{R1_DIR}/seqtab_nochim.rds",
        tsv=f"{R1_DIR}/seqtab_nochim.tsv",
        track=f"{R1_DIR}/track.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        script=f"{workflow.basedir}/scripts/dada2_merge_chimeras.R",
    threads: 16
    resources:
        mem_mb=32000,
        runtime=240,
        tmp_gb=5,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        {params.r_isolation} Rscript {params.script} \
            {threads} {output.rds} {output.tsv} {output.track} {input}
        """


# Same commands as rule qiime2_import
rule qiime2_import_r1:
    input:
        f"{R1_DIR}/seqtab_nochim.tsv",
    output:
        biom=f"{R1_DIR}/feature-table.biom",
        fasta=f"{R1_DIR}/dna-sequences.fasta",
        table=f"{OUT_DIR}/denoise/CFM_ITS1_r1_dada2_table.qza",
        repseqs=f"{OUT_DIR}/denoise/CFM_ITS1_r1_dada2_repseqs.qza",
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


# Same command as rule cluster_its1_otus
rule cluster_its1_r1_otus:
    input:
        table=f"{OUT_DIR}/denoise/CFM_ITS1_r1_dada2_table.qza",
        repseqs=f"{OUT_DIR}/denoise/CFM_ITS1_r1_dada2_repseqs.qza",
    output:
        table=f"{OUT_DIR}/denoise/CFM_ITS1_r1_otu97_table.qza",
        repseqs=f"{OUT_DIR}/denoise/CFM_ITS1_r1_otu97_repseqs.qza",
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


rule qc_its1_r1_vs_merged:
    input:
        r1_table=f"{EXPORT_DIR}/CFM_ITS1_r1_dada2_table/feature-table.biom",
        r1_taxonomy=f"{OUT_DIR}/taxonomy/CFM_ITS1_r1_taxonomy_{MAIN_ITS1}.tsv",
        merged_table=f"{EXPORT_DIR}/CFM_ITS1_dada2_table/feature-table.biom",
        merged_taxonomy=f"{OUT_DIR}/taxonomy/CFM_ITS1_asv_taxonomy_{MAIN_ITS1}.tsv",
    output:
        summary=f"{OUT_DIR}/qc/09_its1_r1_vs_merged.tsv",
        genera=f"{OUT_DIR}/qc/09_its1_r1_vs_merged_genera.tsv",
        figure=f"{OUT_DIR}/qc/09_its1_r1_vs_merged.png",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        script=f"{workflow.basedir}/scripts/qc_its1_r1_vs_merged.py",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        python {params.script} {input.r1_table} {input.r1_taxonomy} \
            {input.merged_table} {input.merged_taxonomy} \
            {output.summary} {output.genera} {output.figure}
        """
