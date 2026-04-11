#!/usr/bin/env python3
"""
Generate QIIME2 PairedEndFastqManifestPhred33V2 manifests for Wemheuer data.

Mirrors the CFM organisation:
    qiime2_wemheuer/import/primer_trimmed/16S/   — trimmed Bac_* reads
    qiime2_wemheuer/import/primer_trimmed/ITS2/  — trimmed Fun_* reads
    qiime2_wemheuer/import/manifest_16S.tsv
    qiime2_wemheuer/import/manifest_ITS2.tsv

Usage (run from project root on the cluster):
    python3 qiime2_wemheuer/scripts/02_make_manifest.py

Inputs:
    data_wemheuer/runinfo.csv              — SRA RunInfo (maps SRR → LibraryName)
    data_wemheuer/raw_data/SRR*_1.fastq.gz — raw paired-end reads from SRA
    qiime2_wemheuer/import/primer_trimmed/16S/SRR*_1.trimmed.fastq.gz   (if trimmed)
    qiime2_wemheuer/import/primer_trimmed/ITS2/SRR*_1.trimmed.fastq.gz  (if trimmed)

Outputs:
    qiime2_wemheuer/import/manifest_16S.tsv
    qiime2_wemheuer/import/manifest_ITS2.tsv
"""

import csv
from pathlib import Path

RUNINFO         = Path("data_wemheuer/runinfo.csv")
RAW_DIR         = Path("data_wemheuer/raw_data")
TRIMMED_16S_DIR = Path("qiime2_wemheuer/import/primer_trimmed/16S")
TRIMMED_ITS2_DIR = Path("qiime2_wemheuer/import/primer_trimmed/ITS2")
OUT_DIR         = Path("qiime2_wemheuer/import")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def find_pair(srr, lib_prefix):
    """Return (r1_abs, r2_abs) for an SRR accession.

    Checks the appropriate primer_trimmed subdirectory first (16S or ITS2),
    then falls back to the raw SRA files.  Returns (None, None) if absent.
    """
    trimmed_dir = TRIMMED_16S_DIR if lib_prefix == "Bac" else TRIMMED_ITS2_DIR

    # trimmed reads produced by cutadapt (Snakefile rule trim_primers_16s/its2)
    r1t = trimmed_dir / f"{srr}_1.trimmed.fastq.gz"
    r2t = trimmed_dir / f"{srr}_2.trimmed.fastq.gz"
    if r1t.exists() and r2t.exists():
        return str(r1t.resolve()), str(r2t.resolve()), "trimmed"

    # raw SRA files (fallback — primers still present)
    r1r = RAW_DIR / f"{srr}_1.fastq.gz"
    r2r = RAW_DIR / f"{srr}_2.fastq.gz"
    if r1r.exists() and r2r.exists():
        return str(r1r.resolve()), str(r2r.resolve()), "raw"

    return None, None, None


def lib_to_sample_id(lib):
    """Map LibraryName → biological sample ID used in QIIME2.

    Fun_T10  → Cocoa_Fungi_10
    Bac_T10  → Cocoa_Bacteria_10
    """
    prefix, _, num = lib.partition("_T")
    n = int(num)
    return f"Cocoa_Fungi_{n:02d}" if prefix == "Fun" else f"Cocoa_Bacteria_{n:02d}"


bac_rows, fun_rows = [], []

with open(RUNINFO) as fh:
    for row in csv.DictReader(fh):
        srr = row["Run"].strip()
        lib = row["LibraryName"].strip()

        if lib.startswith("Bac"):
            lib_prefix = "Bac"
        elif lib.startswith("Fun"):
            lib_prefix = "Fun"
        else:
            print(f"  SKIP (unknown prefix): {srr} ({lib})")
            continue

        r1, r2, src = find_pair(srr, lib_prefix)
        if r1 is None:
            print(f"  MISSING: {srr} ({lib})")
            continue

        sample_id = lib_to_sample_id(lib)
        entry = (sample_id, r1, r2, src)

        if lib_prefix == "Bac":
            bac_rows.append(entry)
        else:
            fun_rows.append(entry)


def write_manifest(path, rows, label):
    with open(path, "w") as fh:
        fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for sample_id, r1, r2, src in sorted(rows):
            fh.write(f"{sample_id}\t{r1}\t{r2}\n")
    print(f"{label}: {len(rows)} samples → {path}")
    # preview first few lines
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i > 3:
                break
            print(f"  {line.rstrip()}")
    # report sources
    sources = {src for *_, src in rows}
    print(f"  Source: {', '.join(sorted(sources))}")


write_manifest(OUT_DIR / "manifest_16S.tsv",  bac_rows, "16S  (Bac)")
write_manifest(OUT_DIR / "manifest_ITS2.tsv", fun_rows, "ITS2 (Fun)")
