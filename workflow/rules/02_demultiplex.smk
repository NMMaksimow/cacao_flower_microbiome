# =============================================================================
# 02. Demultiplexing, per sublibrary and per lane (sequencing run)
# =============================================================================
# Same as v1 step 04 (Stacks process_radtags 2.64, same flags), but each lane on
# its own, and with the v2 barcode files in workflow/resources/internal_tag_mappings,
# which carry the cmt53/cmtPC 16S correction (see the README there).
#
# process_radtags writes one file pair per sample named in the barcode file
# (<sample>_16S.1.fq.gz / .2.fq.gz, plus .rem.* files for reads whose mate was
# discarded) and a log with per-sample counts. The set of samples differs between
# sublibraries, so the whole output folder is the rule's output.
#
# QC: per sample, L1 + L2 decompressed reads, concatenated, must be identical to
# the v1 demultiplexed file (qiime2/import/demux/demultiplexed_sample_files/), and
# cmt53_16S must be identical to the SRA-deposited file.


rule demultiplex:
    input:
        r1=f"{OUT_DIR}/01_adapter_trimmed/{{lane}}/{{sublib}}_R1.fq.gz",
        r2=f"{OUT_DIR}/01_adapter_trimmed/{{lane}}/{{sublib}}_R2.fq.gz",
        barcodes=f"{workflow.basedir}/resources/internal_tag_mappings/{{sublib}}_tags.tsv",
    output:
        directory(f"{OUT_DIR}/02_demux/{{lane}}/{{sublib}}"),
    params:
        base=config["module_base"],
        stacks=config["module_stacks"],
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
        tmp_gb=2,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.stacks}
        mkdir -p {output}
        process_radtags \
            -P \
            -1 {input.r1} \
            -2 {input.r2} \
            -b {input.barcodes} \
            -o {output} \
            --inline_inline \
            --disable-rad-check \
            --retain-header \
            -c -q -r
        """
