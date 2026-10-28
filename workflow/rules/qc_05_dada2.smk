# =============================================================================
# QC 05. Read tracking (main rung) and primer remnants (v1 and all v2 rungs)
# =============================================================================
# Side checks; nothing downstream uses their output.

localrules: qc_read_tracking


rule qc_read_tracking:
    input:
        primer=f"{OUT_DIR}/qc/03_primer_trimming_per_sample.tsv",
        filter_counts=expand(f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/filter_counts.tsv", lane=LANES, allow_missing=True),
        denoise_tracks=expand(
            f"{OUT_DIR}/04_dada2/{{lane}}/{{marker}}/track_{RUNGS[config['main_rung']]['errfun']}.tsv",
            lane=LANES, allow_missing=True,
        ),
        final_track=f"{OUT_DIR}/05_merged/{config['main_rung']}/{{marker}}/track.tsv",
    output:
        f"{OUT_DIR}/qc/05_read_tracking_{{marker}}.tsv",
    params:
        lanes=LANES,
    script:
        "../scripts/qc_read_tracking.py"


rule qc_primer_remnants:
    input:
        v1=expand(
            "qiime2/export/CFM_{marker}_dada2_{part}",
            marker=MARKERS, part=["repseqs/dna-sequences.fasta", "table/feature-table.biom"],
        ),
        v2=expand(f"{OUT_DIR}/05_merged/{{rung}}/{{marker}}/seqtab_nochim.tsv", rung=list(RUNGS), marker=MARKERS),
    output:
        f"{OUT_DIR}/qc/05_primer_remnants.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        script=f"{workflow.basedir}/scripts/qc_primer_remnants.py",
        merged_dir=f"{OUT_DIR}/05_merged",
        rungs=",".join(RUNGS),
        primers=" ".join(f"{m}:{config['primers'][m]['fwd']}:{config['primers'][m]['rev_rc']}" for m in MARKERS),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        python {params.script} {output} {params.merged_dir} qiime2/export {params.rungs} {params.primers}
        """
