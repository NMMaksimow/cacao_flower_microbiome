# =============================================================================
# QC 04. Error-function comparison on mock communities and controls
# =============================================================================
# Side check for the choice of error function: per marker, both runs and both
# functions; nothing downstream uses its output. Expected mock sequences come
# from workflow/resources/mock/ (made by scripts/make_mock_expected.py).

EXPECTED_MOCK = {"16S": "D6305_v3_16S_V4.fasta", "ITS1": "D6305_v3_ITS1.fasta"}


rule qc_error_models:
    input:
        filtered=expand(f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/filtered", lane=LANES, allow_missing=True),
        models=expand(
            f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/errors_{{errfun}}_{{read}}.rds",
            lane=LANES, errfun=config["dada2"]["error_functions"], read=["R1", "R2"], allow_missing=True,
        ),
        expected=lambda wildcards: f"{workflow.basedir}/resources/mock/{EXPECTED_MOCK[wildcards.marker]}",
    output:
        samples=f"{OUT_DIR}/qc/04_error_models_{{marker}}.tsv",
        monotonicity=f"{OUT_DIR}/qc/04_error_models_monotonicity_{{marker}}.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        script=f"{workflow.basedir}/scripts/qc_error_models.R",
        dada2_dir=f"{OUT_DIR}/04_dada2",
        lanes=",".join(LANES),
        errfuns=",".join(config["dada2"]["error_functions"]),
        min_pairs=config["dada2"]["mock_min_pairs"],
    threads: 8
    resources:
        mem_mb=8000,
        runtime=60,
        tmp_gb=2,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        {params.r_isolation} Rscript {params.script} \
            {wildcards.marker} {params.dada2_dir} {params.lanes} {params.errfuns} \
            {input.expected} {params.min_pairs} {threads} \
            {output.samples} {output.monotonicity}
        """
