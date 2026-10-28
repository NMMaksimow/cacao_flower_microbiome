# =============================================================================
# 07. Merge the runs' sequence tables and remove chimeras, per rung and marker
# =============================================================================
# Rung C is the main pipeline; A and B are the ablation rungs (see config
# "ablation_rungs"). Each rung names its lanes and its error function.

RUNGS = config["ablation_rungs"]


def rung_seqtabs(wildcards):
    """Sequence tables of the rung's lanes, denoised with the rung's error function."""
    rung = RUNGS[wildcards.rung]
    return [
        f"{OUT_DIR}/04_dada2/{lane}/{wildcards.marker}/seqtab_{rung['errfun']}.rds"
        for lane in rung["lanes"]
    ]


rule dada2_merge_chimeras:
    input:
        rung_seqtabs,
    output:
        rds=f"{OUT_DIR}/05_merged/{{rung}}/{{marker}}/seqtab_nochim.rds",
        tsv=f"{OUT_DIR}/05_merged/{{rung}}/{{marker}}/seqtab_nochim.tsv",
        track=f"{OUT_DIR}/05_merged/{{rung}}/{{marker}}/track.tsv",
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
