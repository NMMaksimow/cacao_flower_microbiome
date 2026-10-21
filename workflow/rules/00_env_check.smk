# =============================================================================
# 00. Environment check
# =============================================================================
# Loads each software module the pipeline uses and records the tool versions.
# Runs as a short SLURM job, so it also tests the profile (sbatch flags, log
# folder), module loading inside jobs and the R isolation prefix.
# The output is the first part of RUN_INFO.txt.


rule env_check:
    output:
        f"{OUT_DIR}/00_env/tool_versions.txt",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        cutadapt=config["module_cutadapt"],
        stacks=config["module_stacks"],
        r_isolation=config["r_isolation"],
    threads: 1
    resources:
        partition="test",
        mem_mb=2000,
        runtime=5,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base}
        {{
            echo "date: $(date -Iseconds)"
            echo "host: $(hostname)"
            echo "slurm_job_id: ${{SLURM_JOB_ID:-none}}"

            module load {params.cutadapt}
            echo "cutadapt (adapter trimming, {params.cutadapt}): $(cutadapt --version)"
            module unload {params.cutadapt}

            module load {params.stacks}
            echo "stacks ({params.stacks}): $(process_radtags --version 2>&1 | head -n 1)"
            module unload {params.stacks}

            module load {params.qiime2}
            echo "qiime2 ({params.qiime2}): $(qiime --version | head -n 1)"
            echo "cutadapt (primer trimming): $(cutadapt --version)"
            echo "vsearch: $(vsearch --version 2>&1 | head -n 1)"
            {params.r_isolation} Rscript -e 'cat("R:", R.version.string, "\n")' \
                -e 'cat("dada2:", as.character(packageVersion("dada2")), "\n")' \
                -e 'cat("biomformat:", as.character(packageVersion("biomformat")), "\n")'
        }} > {output}
        """
