# =============================================================================
# 04. DADA2 filtering, per lane (sequencing run) and marker
# =============================================================================
# filterAndTrim without truncation (maxEE 2/2, maxN 0, minLen 50, DADA2 defaults
# otherwise). The filtered folder is the output because samples with no surviving
# read pair get no file. Also checks that the reads contain only the quality
# values of the binned error function. Used by the error models (05) and by the
# full DADA2 run.


rule dada2_filter:
    input:
        expand(
            f"{OUT_DIR}/03_primer_trimmed/{{lane}}/{{marker}}/{{sample}}_R{{read}}.fq.gz",
            sample=SAMPLES, read=[1, 2], allow_missing=True,
        ),
    output:
        filtered=directory(f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/filtered"),
        counts=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/filter_counts.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        script=f"{workflow.basedir}/scripts/dada2_filter.R",
        in_dir=lambda wildcards: f"{OUT_DIR}/03_primer_trimmed/{wildcards.lane}/{wildcards.marker}",
        max_ee=" ".join(str(x) for x in config["dada2"]["max_ee"]),
        max_n=config["dada2"]["max_n"],
        min_len=config["dada2"]["min_len"],
        trunc_q=config["dada2"]["trunc_q"],
        qual_bins=",".join(str(q) for q in config["dada2"]["qual_bins"]),
    threads: 16
    resources:
        mem_mb=64000,
        runtime=90,
        tmp_gb=5,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        {params.r_isolation} Rscript {params.script} \
            {params.in_dir} {output.filtered} {output.counts} {threads} \
            {params.max_ee} {params.max_n} {params.min_len} {params.trunc_q} {params.qual_bins}
        """
