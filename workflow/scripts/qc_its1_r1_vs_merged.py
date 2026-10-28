"""QC 09: ITS1 forward reads only (R1) vs merged read pairs, the main ITS1 data.

Pairs from amplicons longer than ~430 nt cannot merge with 2x250 reads; the R1 data
keep them. Both datasets are ASVs classified with the same classifier. Sample types:
mocks (names ending in PC), controls (NC, EC), biosamples (all others).
Reported, as fixed in advance (docs/v2_rerun_plan.md):
  summary  per sample type: reads in R1 and merged, R1/merged per sample (median,
           min, max), fungal share; mocks: share of reads in Saccharomyces and
           Cryptococcus; number of genera lost in the merged data
  genera   per sample type and genus: share of fungal reads in R1 and merged,
           prevalence (samples with the genus) in each; "lost_in_merged" = at least
           1% of fungal R1 reads and absent from the merged data or >= 5-fold lower
  figure   stacked bars of the most abundant genera, R1 vs merged, per sample type
A genus counts unless empty, starting with unidentified/uncultured, containing
Incertae_sedis or ending in "_sp" (as in qc_classifier_comparison.py).

Usage (rule qc_its1_r1_vs_merged, inside the QIIME 2 module for biom/matplotlib):
    python qc_its1_r1_vs_merged.py <r1.biom> <r1_taxonomy.tsv> <merged.biom> \
        <merged_taxonomy.tsv> <summary.tsv> <genera.tsv> <figure.png>
"""
import csv
import sys
from statistics import median

import biom
import matplotlib.pyplot as plt

TYPES = ["mock", "biosample", "control"]
MIN_SHARE = 1.0       # % of fungal R1 reads for a genus to be checked
MIN_FOLD = 5.0        # R1 share / merged share that counts as lost
TOP_PER_TYPE = 6      # most abundant R1 genera per sample type shown in the figure


def sample_type(name):
    if name.endswith("PC"):
        return "mock"
    if name.endswith(("NC", "EC")):
        return "control"
    return "biosample"


def rank(taxon, prefix):
    for part in taxon.split(";"):
        part = part.strip()
        if part.startswith(prefix + "__"):
            name = part[3:]
            if (name and not name.startswith(("unidentified", "uncultured")) and "Incertae_" not in name
                    and not name.endswith("_sp")):
                return name
    return None


def load(table_path, taxonomy_path):
    """Per sample type: {genus or None: reads} for fungal reads, total and fungal reads,
    per-sample totals and per-genus sample sets."""
    table = biom.load_table(table_path)
    with open(taxonomy_path) as handle:
        taxonomy = {row["Feature ID"]: row["Taxon"] for row in csv.DictReader(handle, delimiter="\t")}
    data = {t: {"genera": {}, "total": 0, "fungal": 0, "per_sample": {}, "present": {}} for t in TYPES}
    for values, feature_id, _ in table.iter(axis="observation"):
        taxon = taxonomy[feature_id]
        fungal = rank(taxon, "k") == "Fungi"
        genus = rank(taxon, "g")
        for sample, n in zip(table.ids(), values):
            if n == 0:
                continue
            d = data[sample_type(sample)]
            d["total"] += n
            d["per_sample"][sample] = d["per_sample"].get(sample, 0) + n
            if fungal:
                d["fungal"] += n
                d["genera"][genus] = d["genera"].get(genus, 0) + n
                d["present"].setdefault(genus, set()).add(sample)
    return data


r1_table, r1_taxonomy, merged_table, merged_taxonomy, summary_out, genera_out, figure_out = sys.argv[1:8]
data = {"R1": load(r1_table, r1_taxonomy), "merged": load(merged_table, merged_taxonomy)}


def share(dataset, sample_kind, genus):
    d = data[dataset][sample_kind]
    return 100 * d["genera"].get(genus, 0) / d["fungal"] if d["fungal"] else 0.0


summary, genera_rows = [], []
for t in TYPES:
    r1, merged = data["R1"][t], data["merged"][t]
    ratios = [r1["per_sample"][s] / merged["per_sample"][s] for s in r1["per_sample"] if merged["per_sample"].get(s)]
    lost = 0
    for genus in sorted(set(r1["genera"]) | set(merged["genera"]), key=lambda g: -r1["genera"].get(g, 0)):
        s_r1, s_merged = share("R1", t, genus), share("merged", t, genus)
        is_lost = genus is not None and s_r1 >= MIN_SHARE and (s_merged == 0 or s_r1 / s_merged >= MIN_FOLD)
        lost += is_lost
        genera_rows.append([t, genus or "no genus", f"{s_r1:.3f}", f"{s_merged:.3f}",
                            len(r1["present"].get(genus, ())), len(merged["present"].get(genus, ())), is_lost])
    summary += [
        [t, "samples_R1", len(r1["per_sample"])],
        [t, "reads_R1", int(r1["total"])],
        [t, "reads_merged", int(merged["total"])],
        [t, "R1_over_merged_per_sample_median", f"{median(ratios):.3f}" if ratios else "-"],
        [t, "R1_over_merged_per_sample_min", f"{min(ratios):.3f}" if ratios else "-"],
        [t, "R1_over_merged_per_sample_max", f"{max(ratios):.3f}" if ratios else "-"],
        [t, "fungal_pct_R1", f"{100 * r1['fungal'] / r1['total']:.2f}" if r1["total"] else "-"],
        [t, "fungal_pct_merged", f"{100 * merged['fungal'] / merged['total']:.2f}" if merged["total"] else "-"],
        [t, "genera_lost_in_merged", lost],
    ]
    if t == "mock":
        for genus in ["Saccharomyces", "Cryptococcus"]:
            summary += [[t, f"{genus}_pct_of_fungal_R1", f"{share('R1', t, genus):.2f}"],
                        [t, f"{genus}_pct_of_fungal_merged", f"{share('merged', t, genus):.2f}"]]

with open(summary_out, "w") as out:
    out.write("sample_type\tmetric\tvalue\n")
    for row in summary:
        out.write("\t".join(map(str, row)) + "\n")
with open(genera_out, "w") as out:
    out.write("sample_type\tgenus\tpct_fungal_R1\tpct_fungal_merged\tsamples_R1\tsamples_merged\tlost_in_merged\n")
    for row in genera_rows:
        out.write("\t".join(map(str, row)) + "\n")

# Figure: per sample type, two stacked bars (R1, merged) of genus shares of fungal reads.
# Shown genera: the TOP_PER_TYPE most abundant R1 genera of each sample type (so the
# controls are not hidden behind the biosamples); each keeps one colour in all panels.
# The rest is grouped as "other genera", unassigned as "no genus".
top = []
for t in TYPES:
    named = {g: n for g, n in data["R1"][t]["genera"].items() if g is not None}
    top += [g for g in sorted(named, key=lambda g: -named[g])[:TOP_PER_TYPE] if g not in top]
categories = top + ["other genera", "no genus"]
palette = [c for i, c in enumerate(plt.cm.tab20.colors) if i not in (14, 15)]  # tab20 without its greys
colours = dict(zip(top, palette))                        # one colour per shown genus (max 18)
colours["other genera"] = "#bdbdbd"                      # grey for the rest
colours["no genus"] = "#636363"                          # dark grey for unassigned genus

fig, axes = plt.subplots(1, len(TYPES), figsize=(4 * len(TYPES) + 3, 5), sharey=True)  # panel size in inches
for ax, t in zip(axes, TYPES):
    bottom = [0.0, 0.0]
    for cat in categories:
        heights = []
        for dataset in ["R1", "merged"]:
            if cat == "no genus":
                heights.append(share(dataset, t, None))
            elif cat == "other genera":
                heights.append(sum(share(dataset, t, g) for g in data[dataset][t]["genera"]
                                   if g is not None and g not in top))
            else:
                heights.append(share(dataset, t, cat))
        ax.bar(["R1", "merged"], heights, bottom=bottom, color=colours[cat], label=cat, width=0.6)
        bottom = [b + h for b, h in zip(bottom, heights)]
    ax.set_title(f"{t} (n = {len(data['R1'][t]['per_sample'])})", fontsize=10)
axes[0].set_ylabel("% of fungal reads")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles[::-1], labels[::-1], loc="center right", fontsize=8, frameon=False)  # legend right of panels
fig.suptitle("ITS1: forward reads only (R1) vs merged pairs, genus composition", fontweight="bold")
fig.tight_layout(rect=(0, 0, 0.82, 1))                   # leave the right 18% for the legend
fig.savefig(figure_out, dpi=300)
print("\n".join("\t".join(map(str, row)) for row in summary))
