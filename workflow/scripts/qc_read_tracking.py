"""QC for rule groups 04-07: read pairs per sample at every step (main rung).

Per sample and run: pairs after primer trimming, after filtering, denoised R1 and R2,
merged; then both runs together: merged and non-chimeric, and the share of the
primer-trimmed pairs that ends up in the final table. Samples without a filtered
file in a run (no pair survived filtering) count as 0 from that step on.

Run by Snakemake (script: directive); paths come from the `snakemake` object.
"""
import csv


def read_tsv(path, key="sample"):
    with open(path) as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


marker = snakemake.wildcards.marker
lanes = snakemake.params.lanes
primer = {(row["sample"], row["marker"]): row
          for row in csv.DictReader(open(snakemake.input.primer), delimiter="\t")}
filtered = dict(zip(lanes, (read_tsv(p) for p in snakemake.input.filter_counts)))
denoised = dict(zip(lanes, (read_tsv(p) for p in snakemake.input.denoise_tracks)))
final = read_tsv(snakemake.input.final_track)

samples = sorted(sample for sample, m in primer if m == marker)
columns = ["sample"]
for lane in lanes:
    columns += [f"{lane}_{step}" for step in ("primer_trimmed", "filtered", "denoised_r1", "denoised_r2", "merged")]
columns += ["merged_total", "nonchim", "pct_final_of_trimmed"]

with open(snakemake.output[0], "w") as out:
    out.write("\t".join(columns) + "\n")
    for sample in samples:
        values, trimmed_total = [sample], 0
        for lane in lanes:
            trimmed = int(primer[(sample, marker)][f"{lane}_written"])
            trimmed_total += trimmed
            row = denoised[lane].get(sample, {})
            values += [trimmed, int(filtered[lane][sample]["pairs_out"]),
                       int(float(row.get("denoised_r1", 0))), int(float(row.get("denoised_r2", 0))),
                       int(float(row.get("merged", 0)))]
        row = final.get(sample, {})
        nonchim = int(float(row.get("nonchim", 0)))
        pct = f"{100 * nonchim / trimmed_total:.2f}" if trimmed_total else "NA"
        values += [int(float(row.get("merged", 0))), nonchim, pct]
        out.write("\t".join(map(str, values)) + "\n")
