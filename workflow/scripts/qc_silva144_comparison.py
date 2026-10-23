"""QC 06: compare the four SILVA 144 classifiers on the 16S ASVs (choice criteria).

Criteria fixed in advance (v2 plan): mock recovery; share of reads assigned to genus
on biosamples; no unexpected genera relative to the uniform classifier.
  mock:        taxonomy of the 9 expected mock ASVs (IDs = md5 of the expected
               sequences); genus correct or not; share of mock reads in expected genera
  biosamples:  share of organelle reads (o__Chloroplast, f__Mitochondria) and, among
               bacterial reads (d__Bacteria without organelles, as kept downstream in
               scripts/01a), share with a genus / species assignment; each classifier
               decides itself what is organelle or bacterial (samples whose names do
               not end in PC, EC or NC)
  vs uniform:  ASVs where a weighted classifier names a genus that the uniform
               classifier does not give that ASV (no genus or another genus)
  short ASVs:  class of the ASVs shorter than 200 nt (new in v2)
A rank counts as assigned unless empty or starting with uncultured/Unknown/metagenome
(as in scripts/exploratory/01a_16S_compare_nb_classifiers_CFM.R) or Incertae_Sedis.
SILVA 144 splits some taxa into "X" and "X--other"; both count as X.
The figure shows every metric of the summary table as one bar panel.

Usage (rule qc_silva144_comparison, inside the QIIME 2 module for the biom package):
    python qc_silva144_comparison.py <table.biom> <asv.fasta> <expected.fasta> \
        <min_mock_reads> <summary.tsv> <details.tsv> <figure.png> <name>=<taxonomy.tsv> ...
"""
import csv
import hashlib
import sys
from collections import Counter

import biom
import matplotlib.pyplot as plt

RANKS = ["d", "p", "c", "o", "f", "g", "s"]
UNINFORMATIVE = ("uncultured", "Unknown", "metagenome", "Incertae_Sedis")
# SILVA genus names that are correct for the mock organisms' genus
GENUS_SYNONYMS = {"Escherichia": {"Escherichia-Shigella"}, "Lactobacillus": {"Limosilactobacillus"}}


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            records[name] = records.get(name, "") + line.strip()
    return records


def parse_taxon(taxon):
    """Rank -> name for informative ranks, e.g. {'d': 'Bacteria', ..., 'g': 'Bacillus'}."""
    ranks = {}
    for part in taxon.split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__" and part[0] in RANKS:
            name = part[3:].replace("--other", "")
            if name and not name.startswith(UNINFORMATIVE):
                ranks[part[0]] = name
    return ranks


def is_organelle(ranks):
    return ranks.get("o") == "Chloroplast" or ranks.get("f") == "Mitochondria"


def read_taxonomy(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["Feature ID"]: (parse_taxon(row["Taxon"]), row["Taxon"], float(row["Confidence"]))
                for row in reader}


table_path, asv_fasta, expected_fasta, min_mock_reads, summary_out, details_out, figure_out = sys.argv[1:8]
taxonomies = dict(arg.split("=", 1) for arg in sys.argv[8:])

table = biom.load_table(table_path)
samples = list(table.ids())
sample_sums = dict(zip(samples, table.sum(axis="sample")))
biosamples = [s for s in samples if not s.endswith(("PC", "EC", "NC"))]
mocks = [s for s in samples if s.endswith("PC") and sample_sums[s] >= int(min_mock_reads)]
asv_ids = list(table.ids(axis="observation"))
reads_bio = dict(zip(asv_ids, table.filter(biosamples, inplace=False).sum(axis="observation")))
reads_mock = dict(zip(asv_ids, table.filter(mocks, inplace=False).sum(axis="observation")))
lengths = {i: len(s) for i, s in read_fasta(asv_fasta).items()}

expected = {hashlib.md5(seq.encode()).hexdigest(): name for name, seq in read_fasta(expected_fasta).items()}
expected_genera = {name.split("_")[0] for name in expected.values()}
accepted = set(expected_genera).union(*(GENUS_SYNONYMS.get(g, set()) for g in expected_genera))

tax = {name: read_taxonomy(path) for name, path in taxonomies.items()}
summary, details = [], []
for name, t in tax.items():
    bio_total, mock_total = sum(reads_bio.values()), sum(reads_mock.values())
    genus_ok = 0
    for asv_id, organism in expected.items():
        ranks, taxon, conf = t[asv_id]
        genus = ranks.get("g", "")
        org_genus = organism.split("_")[0]
        ok = genus == org_genus or genus in GENUS_SYNONYMS.get(org_genus, set())
        genus_ok += ok
        details.append([name, "mock_expected_asv", organism.split("|")[0] + "|" + organism.split("|")[1],
                        genus or "-", ranks.get("s", "-"), f"{conf:.3f}", ok])
    mock_in_expected = sum(n for i, n in reads_mock.items() if t[i][0].get("g") in accepted)
    organelle_reads = sum(n for i, n in reads_bio.items() if is_organelle(t[i][0]))
    bacterial = {i: n for i, n in reads_bio.items()
                 if t[i][0].get("d") == "Bacteria" and not is_organelle(t[i][0])}
    bacterial_total = sum(bacterial.values())
    summary += [
        [name, "mock_expected_asvs_genus_correct", f"{genus_ok} of {len(expected)}"],
        [name, "mock_reads_in_expected_genera_pct", f"{100 * mock_in_expected / mock_total:.2f}"],
        [name, "biosample_reads_organelle_pct", f"{100 * organelle_reads / bio_total:.2f}"],
        [name, "bacterial_reads_genus_pct", f"{100 * sum(n for i, n in bacterial.items() if 'g' in t[i][0]) / bacterial_total:.2f}"],
        [name, "bacterial_reads_species_pct", f"{100 * sum(n for i, n in bacterial.items() if 's' in t[i][0]) / bacterial_total:.2f}"],
        [name, "biosample_asvs_genus_pct", f"{100 * sum('g' in t[i][0] for i in asv_ids) / len(asv_ids):.2f}"],
    ]
    if name != "uniform":
        uniform = tax["uniform"]
        unexpected = [i for i in asv_ids
                      if "g" in t[i][0] and t[i][0]["g"] != uniform[i][0].get("g")]
        genus_reads = Counter()
        for i in unexpected:
            genus_reads[t[i][0]["g"]] += reads_bio[i]
        summary += [
            [name, "vs_uniform_asvs_other_genus", len(unexpected)],
            [name, "vs_uniform_biosample_reads_other_genus_pct", f"{100 * sum(reads_bio[i] for i in unexpected) / bio_total:.3f}"],
            [name, "vs_uniform_new_genera_absent_from_uniform",
             len({t[i][0]['g'] for i in unexpected} - {uniform[j][0].get('g') for j in asv_ids})],
        ]
        for genus, n in genus_reads.most_common(10):
            details.append([name, "vs_uniform_top_genera", genus, "-", "-", "-", int(n)])
    short = Counter()
    for i in asv_ids:
        if lengths[i] < 200:
            short[t[i][0].get("c", t[i][0].get("p", t[i][0].get("d", "Unassigned")))] += reads_bio[i]
    for cls, n in short.most_common(8):
        details.append([name, "short_asvs_class", cls, "-", "-", "-", int(n)])

with open(summary_out, "w") as out:
    out.write("classifier\tmetric\tvalue\n")
    for row in summary:
        out.write("\t".join(map(str, row)) + "\n")
with open(details_out, "w") as out:
    out.write("classifier\tcheck\titem\tgenus\tspecies\tconfidence\tvalue\n")
    for row in details:
        out.write("\t".join(map(str, row)) + "\n")

# Figure: one bar panel per summary metric, classifiers on the x axis
metrics = list(dict.fromkeys(row[1] for row in summary))
colours = dict(zip(tax, plt.cm.Set2.colors))           # one colour per classifier (Set2, as in the R plot)
ncols = 3                                              # panels per row
nrows = -(-len(metrics) // ncols)                      # ceiling division
fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))  # panel size in inches
for ax, metric in zip(axes.flat, metrics):
    # "4 of 9" (mock ASVs with correct genus) is plotted as 4
    values = {row[0]: float(str(row[2]).split(" of ")[0]) for row in summary if row[1] == metric}
    bars = ax.bar(list(values), list(values.values()), color=[colours[c] for c in values])
    ax.bar_label(bars, fmt="%.1f", fontsize=8)          # value printed above each bar
    ax.set_title(metric, fontsize=9)                    # metric name as in the summary table
    ax.set_ylabel("%" if metric.endswith("pct") else "count")
    ax.margins(y=0.15)                                  # headroom for the bar labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)  # tilted classifier names
for ax in axes.flat[len(metrics):]:
    ax.axis("off")                                      # hide unused panels
fig.suptitle("QC 06: SILVA 144 classifiers on the 16S ASVs", fontweight="bold")
fig.tight_layout()
fig.savefig(figure_out, dpi=300)
print(f"mocks used: {len(mocks)}, biosamples: {len(biosamples)}, ASVs: {len(asv_ids)}")
