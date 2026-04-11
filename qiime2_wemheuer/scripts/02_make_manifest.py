#!/usr/bin/env python3
"""
Generate QIIME2 PairedEndFastqManifestPhred33V2 manifests for Wemheuer data.
Creates separate manifests for 16S (Bac_*) and ITS1 (Fun_*) libraries.

Usage (run from project root on the cluster):
    python3 qiime2_wemheuer/scripts/02_make_manifest.py

Inputs:
    data_wemheuer/runinfo.csv      — SRA RunInfo (maps SRR → library name)
    data_wemheuer/SRR*_1.fastq.gz — paired-end reads (per sample from SRA)

Outputs:
    qiime2_wemheuer/import/manifest_16S.tsv
    qiime2_wemheuer/import/manifest_ITS1.tsv
"""

import csv
from pathlib import Path

RUNINFO     = Path("data_wemheuer/runinfo.csv")
RAW_DIR     = Path("data_wemheuer")
TRIMMED_DIR = Path("data_wemheuer_trimmed")   # set by trim_primers rules
OUT_DIR     = Path("qiime2_wemheuer/import")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Use trimmed files if available, otherwise fall back to raw SRA files
def find_pair(srr):
    for d, suffix in [(TRIMMED_DIR, ".trimmed.fastq.gz"), (RAW_DIR, ".fastq.gz")]:
        r1 = d.resolve() / f"{srr}_1{suffix}"
        r2 = d.resolve() / f"{srr}_2{suffix}"
        if r1.exists() and r2.exists():
            return str(r1), str(r2), "trimmed" if d == TRIMMED_DIR else "raw"
    return None, None, None

# ── Map SRR → biological sample name ─────────────────────────────────────────
# LibraryName: Fun_T10 → Cocoa_Fungi_10  (matches Table_S1 SeqID)
#              Bac_T10 → Cocoa_Bacteria_10
def lib_to_sample_id(lib):
    prefix, _, num = lib.partition("_T")
    n = int(num)
    return f"Cocoa_Fungi_{n:02d}" if prefix == "Fun" else f"Cocoa_Bacteria_{n:02d}"

bac_rows, fun_rows = [], []

with open(RUNINFO) as fh:
    for row in csv.DictReader(fh):
        srr = row["Run"].strip()
        lib = row["LibraryName"].strip()
        r1, r2, src = find_pair(srr)
        if r1 is None:
            print(f"  MISSING: {srr} ({lib})")
            continue
        sample_id = lib_to_sample_id(lib)
        entry = (sample_id, r1, r2)
        if lib.startswith("Bac"):
            bac_rows.append(entry)
        elif lib.startswith("Fun"):
            fun_rows.append(entry)

print(f"Using {'trimmed' if TRIMMED_DIR.exists() else 'raw'} FASTQs")

def write_manifest(path, rows, label):
    with open(path, "w") as fh:
        fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for sample_id, r1, r2 in sorted(rows):
            fh.write(f"{sample_id}\t{r1}\t{r2}\n")
    print(f"{label}: {len(rows)} samples → {path}")
    # preview
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i > 3: break
            print(f"  {line.rstrip()}")

write_manifest(OUT_DIR / "manifest_16S.tsv",  bac_rows, "16S (Bac)")
write_manifest(OUT_DIR / "manifest_ITS1.tsv", fun_rows, "ITS1 (Fun)")
