# =============================================================================
# 06. DADA2 denoising and read merging, per lane, marker and error model
# =============================================================================
# lane = L1 or L2 (one sequencing run each), or "pooled" (both runs' filtered
# reads concatenated per sample, as v1 did; ablation rung A only). The error
# model comes from rule dada2_learn_errors (05), which works for "pooled" too.


rule dada2_pool_runs:
    input:
        l1=f"{OUT_DIR}/04_dada2/L1/{{marker}}/filtered",
        l2=f"{OUT_DIR}/04_dada2/L2/{{marker}}/filtered",
    output:
        directory(f"{OUT_DIR}/04_dada2/pooled/{{marker}}/filtered"),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output}
        # A sample may have no file in one run if no read pair survived filtering
        for name in $(ls {input.l1} {input.l2} | grep "\.fq\.gz$" | sort -u); do
            for run_dir in {input.l1} {input.l2}; do
                if [ -e "$run_dir/$name" ]; then cat "$run_dir/$name"; fi
            done > {output}/$name
        done
        """


rule dada2_denoise:
    input:
        filtered=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/filtered",
        err_r1=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_R1.rds",
        err_r2=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_R2.rds",
    output:
        seqtab=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/seqtab_{{errfun}}.rds",
        track=f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/track_{{errfun}}.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        script=f"{workflow.basedir}/scripts/dada2_denoise.R",
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
            {input.filtered} {input.err_r1} {input.err_r2} {threads} \
            {output.seqtab} {output.track}
        """
