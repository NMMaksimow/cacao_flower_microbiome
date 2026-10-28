# =============================================================================
# 03. Primer trimming, per lane (sequencing run) and marker
# =============================================================================
# Replaces v1 step 05, whose 3' options were given the primers instead of their
# reverse complements. One linked adapter per read (cutadapt recipe for reads
# that can be longer than the amplicon):
#   R1: ^FWD...RC(REV)   R2: ^REV...RC(FWD)
# The anchored 5' primer is required; the 3' part is optional and, when found,
# is removed with everything after it (inline tag, Nextera adapter). Pairs where
# R1 or R2 does not start with its primer are discarded. The same command serves
# 16S and ITS1: long ITS1 inserts have no read-through and keep their full length.
#
# One job per lane x marker loops over all samples (each takes seconds); the
# output files are still listed per sample. The .rem files of demultiplexing
# (reads whose mate was removed) are not used, as in v1.


def sample_sublibs(wildcards):
    """'sample:sublibrary' pairs for the shell loop of one marker."""
    return " ".join(f"{s}:{SUBLIB_OF[(s, wildcards.marker)]}" for s in SAMPLES)


rule primer_trimming:
    input:
        expand(f"{OUT_DIR}/02_demux/{{lane}}/{{sublib}}", sublib=SUBLIBS, allow_missing=True),
    output:
        r1=expand(f"{OUT_DIR}/03_primer_trimmed/{{lane}}/{{marker}}/{{sample}}_R1.fq.gz", sample=SAMPLES, allow_missing=True),
        r2=expand(f"{OUT_DIR}/03_primer_trimmed/{{lane}}/{{marker}}/{{sample}}_R2.fq.gz", sample=SAMPLES, allow_missing=True),
        json=expand(f"{OUT_DIR}/03_primer_trimmed/{{lane}}/{{marker}}/{{sample}}.cutadapt.json", sample=SAMPLES, allow_missing=True),
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        samples=sample_sublibs,
        demux_dir=lambda wildcards: f"{OUT_DIR}/02_demux/{wildcards.lane}",
        out_dir=lambda wildcards: f"{OUT_DIR}/03_primer_trimmed/{wildcards.lane}/{wildcards.marker}",
        fwd=lambda wildcards: config["primers"][wildcards.marker]["fwd"],
        rev=lambda wildcards: config["primers"][wildcards.marker]["rev"],
        fwd_rc=lambda wildcards: config["primers"][wildcards.marker]["fwd_rc"],
        rev_rc=lambda wildcards: config["primers"][wildcards.marker]["rev_rc"],
        min_length=config["primer_min_length"],
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        tmp_gb=1,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        for pair in {params.samples}; do
            sample=${{pair%%:*}}
            sublib=${{pair#*:}}
            echo "== $sample ({wildcards.marker}, {wildcards.lane}, $sublib)"
            cutadapt \
                -a "^{params.fwd}...{params.rev_rc}" \
                -A "^{params.rev}...{params.fwd_rc}" \
                --match-read-wildcards \
                --discard-untrimmed \
                --pair-filter=any \
                --minimum-length {params.min_length} \
                --report=minimal \
                --json "{params.out_dir}/$sample.cutadapt.json" \
                -o "{params.out_dir}/${{sample}}_R1.fq.gz" \
                -p "{params.out_dir}/${{sample}}_R2.fq.gz" \
                "{params.demux_dir}/$sublib/${{sample}}_{wildcards.marker}.1.fq.gz" \
                "{params.demux_dir}/$sublib/${{sample}}_{wildcards.marker}.2.fq.gz"
        done
        """
