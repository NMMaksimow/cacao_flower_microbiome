"""QC: primer remnants and length distribution of the final ASVs, v1 vs v2 rungs.

For every ASV of v1 (qiime2/export, read only) and of the v2 rungs A, B, C:
  fwd_start    forward primer within the first 15 nt (+ primer length), i.e. at
               the 5' end, possibly behind an inline-tag remnant;
  rcrev_any    reverse complement of the reverse primer anywhere (3' read-through);
  rcrev_end    a partial reverse-complement reverse primer (10 nt or more) at the
               3' end, with at most 1 mismatch per 10 nt.
Full primers may differ by up to 2 mismatches (cutadapt's tolerance). Results are
given as numbers of ASVs and as shares of reads, plus read-weighted ASV lengths.

Usage (rule qc_primer_remnants, inside the QIIME 2 module for the biom package):
    python qc_primer_remnants.py <out.tsv> <merged_dir> <v1_export_dir> <rungs> \
        <marker:FWD:RC_REV> [<marker:FWD:RC_REV> ...]
"""
import csv
import re
import sys

import biom

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
         "V": "ACG", "N": "ACGT"}
MAX_MISMATCHES = 2
START_WINDOW = 15


def mismatches(primer, site):
    return sum(base not in IUPAC[p] for p, base in zip(primer, site))


def primer_sites(seq, primer):
    """Start positions where the primer matches with <= MAX_MISMATCHES.

    With at most 2 mismatches one third of the primer matches exactly, so exact
    matches of the thirds are used as seeds before the full comparison.
    """
    third = len(primer) // 3
    seeds = [(0, primer[:third]), (third, primer[third:2 * third]), (2 * third, primer[2 * third:])]
    starts = set()
    for offset, part in seeds:
        for match in re.finditer("".join(f"[{IUPAC[b]}]" for b in part), seq):
            start = match.start() - offset
            if 0 <= start <= len(seq) - len(primer) and mismatches(primer, seq[start:start + len(primer)]) <= MAX_MISMATCHES:
                starts.add(start)
    return starts


def remnants(seq, fwd, rc_rev):
    return {
        "fwd_start": any(start <= START_WINDOW for start in primer_sites(seq, fwd)),
        "rcrev_any": bool(primer_sites(seq, rc_rev)),
        "rcrev_end": any(mismatches(rc_rev[:k], seq[-k:]) <= k // 10 for k in range(10, len(rc_rev))),
    }


def weighted_quantile(lengths, weights, q):
    total, cumulative = sum(weights), 0
    for length, weight in sorted(zip(lengths, weights)):
        cumulative += weight
        if cumulative >= q * total:
            return length


def v1_asvs(export_dir, marker):
    """(sequence, reads) of the v1 ASV table of one marker."""
    seqs, name = {}, None
    for line in open(f"{export_dir}/CFM_{marker}_dada2_repseqs/dna-sequences.fasta"):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            seqs[name] = line.strip()
    table = biom.load_table(f"{export_dir}/CFM_{marker}_dada2_table/feature-table.biom")
    sums = table.sum(axis="observation")
    return [(seqs[i], int(n)) for i, n in zip(table.ids(axis="observation"), sums)]


def v2_asvs(merged_dir, rung, marker):
    """(sequence, reads) of one v2 rung, from its sequence-table TSV."""
    with open(f"{merged_dir}/{rung}/{marker}/seqtab_nochim.tsv") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        return [(row[0], sum(int(float(x)) for x in row[1:])) for row in reader]


out_tsv, merged_dir, v1_dir, rungs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
primers = {m: (f, r) for m, f, r in (arg.split(":") for arg in sys.argv[5:])}

checks = ["fwd_start", "rcrev_any", "rcrev_end", "any_remnant"]
header = ["source", "marker", "asvs", "reads"]
header += [f"asvs_{c}" for c in checks] + [f"pct_reads_{c}" for c in checks]
header += ["len_min", "len_p01", "len_median", "len_max"]

with open(out_tsv, "w") as out:
    out.write("\t".join(header) + "\n")
    for marker, (fwd, rc_rev) in primers.items():
        sources = [("v1", v1_asvs(v1_dir, marker))] + [(f"v2_{r}", v2_asvs(merged_dir, r, marker)) for r in rungs]
        for source, asvs in sources:
            hits = [remnants(seq, fwd, rc_rev) for seq, _ in asvs]
            for h in hits:
                h["any_remnant"] = h["fwd_start"] or h["rcrev_any"] or h["rcrev_end"]
            reads = [n for _, n in asvs]
            lengths = [len(seq) for seq, _ in asvs]
            total = sum(reads)
            row = [source, marker, len(asvs), total]
            row += [sum(h[c] for h in hits) for c in checks]
            row += [f"{100 * sum(n for h, n in zip(hits, reads) if h[c]) / total:.3f}" for c in checks]
            row += [min(lengths), weighted_quantile(lengths, reads, 0.01),
                    weighted_quantile(lengths, reads, 0.5), max(lengths)]
            out.write("\t".join(map(str, row)) + "\n")
            print(source, marker, "done")
