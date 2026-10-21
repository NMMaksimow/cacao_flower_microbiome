# =============================================================================
# QC 03. Primer trimming summary (v2 per lane and marker, v1 for comparison)
# =============================================================================
# Side check: reads the cutadapt JSON reports (and the v1 log) and writes two
# small tables. It takes seconds, so it runs on the controller instead of as a
# SLURM job ("localrules").

localrules: qc_primer_trimming


rule qc_primer_trimming:
    input:
        expand(
            f"{OUT_DIR}/03_primer_trimmed/{{lane}}/{{marker}}/{{sample}}.cutadapt.json",
            lane=LANES,
            marker=MARKERS,
            sample=SAMPLES,
        ),
    output:
        summary=f"{OUT_DIR}/qc/03_primer_trimming.tsv",
        samples=f"{OUT_DIR}/qc/03_primer_trimming_per_sample.tsv",
    params:
        v1_log="logs/05_primer_trimming.out",
        primers=config["primers"],
    script:
        "../scripts/qc_primer_trimming.py"
