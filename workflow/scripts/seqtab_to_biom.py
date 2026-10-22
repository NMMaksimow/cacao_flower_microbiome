"""Convert a DADA2 sequence table (TSV) into a BIOM table and a FASTA file.

Input: TSV with one row per ASV (column "sequence", then one column per sample).
ASV IDs are the md5 hash of the sequence, the QIIME 2 convention (v1 used it too),
so the same sequence gets the same ID in every run and table.

Usage (rule qiime2_import):
    python seqtab_to_biom.py <seqtab_nochim.tsv> <feature-table.biom> <dna-sequences.fasta>
"""
import csv
import hashlib
import sys

import biom
import numpy as np
from biom.util import biom_open

tsv_in, biom_out, fasta_out = sys.argv[1:4]

with open(tsv_in) as handle:
    reader = csv.reader(handle, delimiter="\t")
    samples = next(reader)[1:]
    rows = [(row[0], [int(float(x)) for x in row[1:]]) for row in reader]

ids = [hashlib.md5(seq.encode()).hexdigest() for seq, _ in rows]
table = biom.Table(np.array([counts for _, counts in rows]), observation_ids=ids, sample_ids=samples)
with biom_open(biom_out, "w") as handle:
    table.to_hdf5(handle, "CFM v2 DADA2")

with open(fasta_out, "w") as out:
    for asv_id, (seq, _) in zip(ids, rows):
        out.write(f">{asv_id}\n{seq}\n")

print(f"{len(ids)} ASVs x {len(samples)} samples, {int(table.sum())} reads")
