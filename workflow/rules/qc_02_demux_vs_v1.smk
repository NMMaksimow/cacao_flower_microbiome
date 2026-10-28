# =============================================================================
# QC 02. v2 demultiplexing vs v1 and SRA
# =============================================================================
# Side check, not a pipeline step: nothing downstream uses its output.
# For every v1 demultiplexed file (2 x 1,344), the decompressed content must
# equal v2 L1 followed by v2 L2; cmt53_16S must also equal the SRA copy.
# Output: one row per reference file, with md5 values and identical = True/False.


rule qc_demux_vs_v1:
    input:
        expand(f"{OUT_DIR}/02_demux/{{lane}}/{{sublib}}", lane=LANES, sublib=SUBLIBS),
    output:
        f"{OUT_DIR}/qc/02_demux_vs_v1.tsv",
    params:
        demux_dir=f"{OUT_DIR}/02_demux",
        lanes=LANES,
        v1_dir="qiime2/import/demux/demultiplexed_sample_files",
        sra_dir="sra_upload/16S",
    threads: 16
    resources:
        mem_mb=8000,
        runtime=60,
        tmp_gb=1,
    script:
        "../scripts/qc_demux_vs_v1.py"
