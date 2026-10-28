# =============================================================================
# 05. DADA2 error models, per lane (sequencing run), marker and error function
# =============================================================================
# Each run gets its own error model (the two lanes are different flowcells and
# instruments). Two error functions are learned for comparison:
#   loess  = DADA2 default loessErrfun
#   binned = makeBinnedQualErrfun at the NovaSeq quality bins
# learnErrors uses randomly chosen samples (fixed seed) up to learn_nbases bases.


rule dada2_learn_errors:
    input:
        f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/filtered",
    output:
        rds_r1=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_R1.rds",
        rds_r2=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_R2.rds",
        png_r1=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_R1.png",
        png_r2=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_R2.png",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        script=f"{workflow.basedir}/scripts/dada2_learn_errors.R",
        qual_bins=",".join(str(q) for q in config["dada2"]["qual_bins"]),
        nbases=config["dada2"]["learn_nbases"],
        randomize="TRUE" if config["dada2"]["learn_randomize"] else "FALSE",
        seed=config["dada2"]["seed"],
    threads: 16
    resources:
        mem_mb=8000,
        runtime=120,
        tmp_gb=5,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        {params.r_isolation} Rscript {params.script} \
            {input} {wildcards.errfun} {params.qual_bins} {params.nbases} \
            {params.randomize} {params.seed} {threads} \
            {output.rds_r1} {output.rds_r2} {output.png_r1} {output.png_r2} \
            "{wildcards.lane} | {wildcards.marker}"
        """
