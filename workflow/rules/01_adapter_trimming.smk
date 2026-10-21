# =============================================================================
# 01. Adapter trimming, per sublibrary and per lane (sequencing run)
# =============================================================================
# Same as v1 step 02 (cutadapt 4.9, same Nextera adapters, default settings),
# but each lane is trimmed on its own instead of L1+L2 concatenated. cutadapt
# treats every read independently, so the trimmed reads equal those of v1.
# QC: per sublibrary, L1 + L2 counts in the JSON reports must equal the counts
# in the v1 log (logs/02_adapter_trimming.out).


def raw_fastq(wildcards, read):
    """Path of one raw FASTQ file, built from the sublibrary, lane and read."""
    flowcell = config["flowcells"][wildcards.lane]
    return (
        f"{RAW_DIR}/{wildcards.sublib}/"
        f"{wildcards.sublib}_{config['library_id']}_{flowcell}_{wildcards.lane}_{read}.fq.gz"
    )


rule adapter_trimming:
    input:
        r1=lambda wildcards: raw_fastq(wildcards, 1),
        r2=lambda wildcards: raw_fastq(wildcards, 2),
    output:
        r1=f"{OUT_DIR}/01_adapter_trimmed/{{lane}}/{{sublib}}_R1.fq.gz",
        r2=f"{OUT_DIR}/01_adapter_trimmed/{{lane}}/{{sublib}}_R2.fq.gz",
        json=f"{OUT_DIR}/01_adapter_trimmed/{{lane}}/{{sublib}}.cutadapt.json",
    params:
        base=config["module_base"],
        cutadapt=config["module_cutadapt"],
        adapter_r1=config["adapter_r1"],
        adapter_r2=config["adapter_r2"],
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.cutadapt}
        cutadapt \
            -a {params.adapter_r1} \
            -A {params.adapter_r2} \
            --json {output.json} \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2}
        """
