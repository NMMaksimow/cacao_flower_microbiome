"""QC 06 / 07: compare taxonomy classifiers on the same features (choice criteria).

16S (QC 06): the four SILVA 144 classifiers on the ASVs.
ITS1 (QC 07): the two UNITE classifiers on the 97% OTUs.
Criteria fixed in advance (v2 plan): mock recovery; share of target reads assigned to
genus on biosamples; few genera that contradict the reference classifier.
  mock:        taxonomy of the expected mock sequences that are features (IDs = md5 of
               the expected sequences; for OTUs only centroids can match); genus
               correct or not; expected genera found; share of mock reads in them
  biosamples:  share of target reads, and among them the share with a genus / species
               assignment. Target = what is kept downstream: 16S d__Bacteria without
               organelles (scripts/01a), ITS1 k__Fungi (scripts/01b); each classifier
               decides itself what is target (samples whose names do not end in PC,
               EC or NC). 16S also: share of organelle reads
  vs reference: features where a classifier names a genus that the reference (the
               first classifier given: 16S uniform, ITS1 s_all = v1's) does not give
               that feature (no genus or another genus)
  short ASVs:  16S only, class of the ASVs shorter than 200 nt (new in v2)
A rank counts as assigned unless empty, starting with uncultured/Unknown/metagenome/
unidentified, containing Incertae_Sedis/Incertae_sedis anywhere (SILVA "Incertae_Sedis";
UNITE placeholders such as "Pleosporaceae_gen_Incertae_sedis") or ending in "_sp"
(UNITE species placeholders such as "Fungi_sp").
SILVA 144 splits some taxa into "X" and "X--other"; both count as X.
The figure shows every metric of the summary table as one bar panel.

Usage (rules qc_silva144_comparison, qc_unite_comparison; inside the QIIME 2 module
for the biom package):
    python qc_classifier_comparison.py <marker> <table.biom> <features.fasta> \
        <expected.fasta> <min_mock_reads> <summary.tsv> <details.tsv> <figure.png> \
        <reference_name>=<taxonomy.tsv> <name>=<taxonomy.tsv> ...
"""
import csv
import hashlib
import sys
from collections import Counter

import biom
import matplotlib.pyplot as plt

RANKS = ["d", "k", "p", "c", "o", "f", "g", "s"]
UNINFORMATIVE = ("uncultured", "Unknown", "metagenome", "unidentified", "Incertae_Sedis", "Incertae_sedis")
# Target reads per marker: (rank, name) kept downstream
TARGET = {"16S": ("d", "Bacteria"), "ITS1": ("k", "Fungi")}
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
            # UNITE placeholders: "Pleosporaceae_gen_Incertae_sedis", "Fungi_sp"
            if (name and not name.startswith(UNINFORMATIVE) and "Incertae_" not in name
                    and not name.endswith("_sp")):
                ranks[part[0]] = name
    return ranks


def is_organelle(ranks):
    return ranks.get("o") == "Chloroplast" or ranks.get("f") == "Mitochondria"


def is_target(ranks, marker):
    rank, name = TARGET[marker]
    return ranks.get(rank) == name and not is_organelle(ranks)


def read_taxonomy(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["Feature ID"]: (parse_taxon(row["Taxon"]), row["Taxon"], float(row["Confidence"]))
                for row in reader}


marker, table_path, features_fasta, expected_fasta, min_mock_reads, summary_out, details_out, figure_out = sys.argv[1:9]
taxonomies = dict(arg.split("=", 1) for arg in sys.argv[9:])
reference = next(iter(taxonomies))

table = biom.load_table(table_path)
samples = list(table.ids())
sample_sums = dict(zip(samples, table.sum(axis="sample")))
biosamples = [s for s in samples if not s.endswith(("PC", "EC", "NC"))]
mocks = [s for s in samples if s.endswith("PC") and sample_sums[s] >= int(min_mock_reads)]
feature_ids = list(table.ids(axis="observation"))
reads_bio = dict(zip(feature_ids, table.filter(biosamples, inplace=False).sum(axis="observation")))
reads_mock = dict(zip(feature_ids, table.filter(mocks, inplace=False).sum(axis="observation")))
lengths = {i: len(s) for i, s in read_fasta(features_fasta).items()}

expected = {hashlib.md5(seq.encode()).hexdigest(): name for name, seq in read_fasta(expected_fasta).items()}
expected_genera = {name.split("_")[0] for name in expected.values()}
accepted = set(expected_genera).union(*(GENUS_SYNONYMS.get(g, set()) for g in expected_genera))
expected_present = [i for i in expected if i in reads_mock]

tax = {name: read_taxonomy(path) for name, path in taxonomies.items()}
summary, details = [], []
for name, t in tax.items():
    bio_total, mock_total = sum(reads_bio.values()), sum(reads_mock.values())
    genus_ok = 0
    for feature_id in expected_present:
        ranks, taxon, conf = t[feature_id]
        organism = expected[feature_id]
        genus = ranks.get("g", "")
        org_genus = organism.split("_")[0]
        ok = genus == org_genus or genus in GENUS_SYNONYMS.get(org_genus, set())
        genus_ok += ok
        details.append([name, "mock_expected_feature", organism.split("|")[0] + "|" + organism.split("|")[1],
                        genus or "-", ranks.get("s", "-"), f"{conf:.3f}", ok])
    mock_genera = {t[i][0].get("g") for i, n in reads_mock.items() if n > 0}
    genera_found = sum(bool(({g} | GENUS_SYNONYMS.get(g, set())) & mock_genera) for g in expected_genera)
    mock_in_expected = sum(n for i, n in reads_mock.items() if t[i][0].get("g") in accepted)
    target = {i: n for i, n in reads_bio.items() if is_target(t[i][0], marker)}
    target_total = sum(target.values())
    summary += [
        [name, "mock_expected_features_genus_correct", f"{genus_ok} of {len(expected_present)}"],
        [name, "mock_expected_genera_found", f"{genera_found} of {len(expected_genera)}"],
        [name, "mock_reads_in_expected_genera_pct", f"{100 * mock_in_expected / mock_total:.2f}"],
        [name, "biosample_reads_target_pct", f"{100 * target_total / bio_total:.2f}"],
        [name, "target_reads_genus_pct", f"{100 * sum(n for i, n in target.items() if 'g' in t[i][0]) / target_total:.2f}"],
        [name, "target_reads_species_pct", f"{100 * sum(n for i, n in target.items() if 's' in t[i][0]) / target_total:.2f}"],
        [name, "biosample_features_genus_pct", f"{100 * sum('g' in t[i][0] for i in feature_ids) / len(feature_ids):.2f}"],
    ]
    if marker == "16S":
        organelle_reads = sum(n for i, n in reads_bio.items() if is_organelle(t[i][0]))
        summary.append([name, "biosample_reads_organelle_pct", f"{100 * organelle_reads / bio_total:.2f}"])
    if name != reference:
        ref = tax[reference]
        unexpected = [i for i in feature_ids
                      if "g" in t[i][0] and t[i][0]["g"] != ref[i][0].get("g")]
        contradicting = [i for i in unexpected if "g" in ref[i][0]]
        genus_reads = Counter()
        for i in unexpected:
            genus_reads[t[i][0]["g"]] += reads_bio[i]
        summary += [
            [name, f"vs_{reference}_features_other_genus", len(unexpected)],
            [name, f"vs_{reference}_biosample_reads_other_genus_pct", f"{100 * sum(reads_bio[i] for i in unexpected) / bio_total:.3f}"],
            [name, f"vs_{reference}_target_reads_contradicting_genus_pct",
             f"{100 * sum(target.get(i, 0) for i in contradicting) / target_total:.3f}"],
            [name, f"vs_{reference}_new_genera_absent_from_{reference}",
             len({t[i][0]['g'] for i in unexpected} - {ref[j][0].get('g') for j in feature_ids})],
        ]
        for genus, n in genus_reads.most_common(10):
            details.append([name, f"vs_{reference}_top_genera", genus, "-", "-", "-", int(n)])
    if marker == "16S":
        short = Counter()
        for i in feature_ids:
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
    # "4 of 9" (mock features with correct genus) is plotted as 4
    values = {row[0]: float(str(row[2]).split(" of ")[0]) for row in summary if row[1] == metric}
    bars = ax.bar(list(values), list(values.values()), color=[colours[c] for c in values])
    ax.bar_label(bars, fmt="%.1f", fontsize=8)          # value printed above each bar
    ax.set_title(metric, fontsize=9)                    # metric name as in the summary table
    ax.set_ylabel("%" if metric.endswith("pct") else "count")
    ax.margins(y=0.15)                                  # headroom for the bar labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)  # tilted classifier names
for ax in axes.flat[len(metrics):]:
    ax.axis("off")                                      # hide unused panels
fig.suptitle(f"{marker}: classifier comparison on {len(feature_ids)} features", fontweight="bold")
fig.tight_layout()
fig.savefig(figure_out, dpi=300)
print(f"mocks used: {len(mocks)}, biosamples: {len(biosamples)}, features: {len(feature_ids)}")
