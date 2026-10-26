"""QC 10: v1 vs v2 comparison table, with the ablation rungs (Phase 2.8).

Datasets: v1 (qiime2/, q2-dada2 2024.10 on lanes merged before DADA2, no 3' primer
trimming) and the v2 rungs A (runs pooled, loess error model), B (per run, loess) and
C (per run, binned = final v2). For top taxa also "v1_reclassified": v1 sequences
with the v2 main classifiers (rule file 16), which separates the classifier effect
(v1 -> v1_reclassified) from the processing effect (v1_reclassified -> C).
A and B have no taxonomy and no OTUs ("-").

Per marker (rows = metrics, columns = datasets):
  retention    read pairs summed over all samples / over biosamples after each step
               (demultiplexed, primer-trimmed, filtered, denoised, merged,
               non-chimeric), as % of demultiplexed; median per-biosample % kept.
               v1 did not discard reads at primer trimming (primer-trimmed = input);
               "denoised" counts forward reads, as q2-dada2 does
  features     ASVs, ASVs present in biosamples, median ASV length, ITS1 97% OTUs
  remnants     % of reads with a primer remnant, 1st-percentile ASV length (QC 05)
  mock         mocks with >= min_mock_reads read pairs: expected sequences found,
               % of mock reads in them, share per organism, other ASVs with >= 0.1%
               of mock reads. "exact" = identical sequence; "contained" = the
               expected sequence lies inside the ASV (v1 ASVs may carry 3' primer
               and tag remnants)
  top_genera   the 10 most abundant genera of each of v1, v1_reclassified and C in
               biosample target reads (16S: Bacteria without chloroplast and
               mitochondria, ASVs; ITS1: Fungi, 97% OTUs), % of target reads
Figure: rows 16S / ITS1; columns: retention per step, mock composition, top genera.

Usage (rule qc_v1_vs_v2, inside the QIIME 2 module for biom/matplotlib; working
directory = project root, paths follow the pipeline layout):
    python qc_v1_vs_v2.py <v2_dir> <v1_dir> <mock_dir> <min_mock_reads> <rungs> \
        <main_16S> <main_ITS1> <out_16S.tsv> <out_ITS1.tsv> <figure.png>
    <rungs>  e.g. "A=pooled:loess,B=L1+L2:loess,C=L1+L2:binned" (config ablation_rungs)
"""
import csv
import hashlib
import sys
import zipfile
from collections import Counter
from statistics import median

import biom
import matplotlib.pyplot as plt

(v2_dir, v1_dir, mock_dir, min_mock, rung_spec, main_16S, main_ITS1,
 out_16S, out_ITS1, figure_out) = sys.argv[1:11]
min_mock = int(min_mock)
RUNGS = {name: (lanes.split("+"), errfun)
         for name, spec in (item.split("=") for item in rung_spec.split(","))
         for lanes, errfun in [spec.split(":")]}
DATASETS = ["v1"] + list(RUNGS)
STEPS = ["demultiplexed", "primer_trimmed", "filtered", "denoised", "merged", "nonchim"]
MOCK_FILES = {"16S": "D6305_v3_16S_V4.fasta", "ITS1": "D6305_v3_ITS1.fasta"}
UNINFORMATIVE = ("uncultured", "Unknown", "metagenome", "unidentified")


def is_biosample(sample):
    return not sample.endswith(("PC", "EC", "NC"))


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            records[name] = records.get(name, "") + line.strip()
    return records


def genus_and_target(taxon, marker):
    """(genus or None, is target read) from a QIIME 2 taxonomy string."""
    ranks = {}
    for part in taxon.split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__":
            name = part[3:].replace("--other", "")
            if (name and not name.startswith(UNINFORMATIVE) and "Incertae_" not in name
                    and not name.endswith("_sp")):
                ranks[part[0]] = name
    if marker == "16S":
        target = (ranks.get("d") == "Bacteria" and ranks.get("o") != "Chloroplast"
                  and ranks.get("f") != "Mitochondria")
    else:
        target = ranks.get("k") == "Fungi"
    return ranks.get("g"), target


def read_taxonomy(path, marker):
    return {row["Feature ID"]: genus_and_target(row["Taxon"], marker) for row in read_tsv(path)}


# ── Features: (id, sequence, biosample reads, {mock sample: reads}) per dataset ──
def features_from_seqtab(path):
    """Rows of a DADA2 seqtab_nochim.tsv (sequence + one column per sample)."""
    out = []
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        samples = next(reader)[1:]
        bio = [i for i, s in enumerate(samples) if is_biosample(s)]
        mock = [i for i, s in enumerate(samples) if s.endswith("PC")]
        for row in reader:
            counts = [int(float(x)) for x in row[1:]]
            out.append((hashlib.md5(row[0].encode()).hexdigest(), row[0],
                        sum(counts[i] for i in bio),
                        {samples[i]: counts[i] for i in mock if counts[i]}))
    return out


def features_from_biom(table_path, fasta_path=None):
    table = biom.load_table(table_path)
    samples = list(table.ids())
    seqs = read_fasta(fasta_path) if fasta_path else {}
    out = []
    for values, feature_id, _ in table.iter(axis="observation"):
        counts = dict(zip(samples, values))
        out.append((feature_id, seqs.get(feature_id, ""),
                    int(sum(n for s, n in counts.items() if is_biosample(s))),
                    {s: int(n) for s, n in counts.items() if s.endswith("PC") and n}))
    return out


# ── Read retention ──
def v1_stats(marker):
    """Per-sample q2-dada2 statistics, read from the artifact without unpacking."""
    with zipfile.ZipFile(f"{v1_dir}/denoise/CFM_{marker}_dada2_stats.qza") as archive:
        name = next(n for n in archive.namelist() if n.endswith("/data/stats.tsv"))
        lines = [l for l in archive.read(name).decode().splitlines() if not l.startswith("#q2")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    per_step = {step: {} for step in STEPS}
    for row in rows:
        s = row["sample-id"]
        for step, column in [("demultiplexed", "input"), ("primer_trimmed", "input"),
                             ("filtered", "filtered"), ("denoised", "denoised"),
                             ("merged", "merged"), ("nonchim", "non-chimeric")]:
            per_step[step][s] = int(float(row[column]))
    return per_step


def v2_steps(marker, rung):
    lanes, errfun = RUNGS[rung]
    per_step = {step: Counter() for step in STEPS}
    for row in read_tsv(f"{v2_dir}/qc/03_primer_trimming_per_sample.tsv"):
        if row["marker"] == marker:
            per_step["demultiplexed"][row["sample"]] = int(row["L1_in"] or 0) + int(row["L2_in"] or 0)
            per_step["primer_trimmed"][row["sample"]] = int(row["L1_written"] or 0) + int(row["L2_written"] or 0)
    for lane in ["L1", "L2"]:
        for row in read_tsv(f"{v2_dir}/04_dada2/{lane}/{marker}/filter_counts.tsv"):
            per_step["filtered"][row["sample"]] += int(row["pairs_out"])
    for lane in lanes:
        for row in read_tsv(f"{v2_dir}/04_dada2/{lane}/{marker}/track_{errfun}.tsv"):
            per_step["denoised"][row["sample"]] += int(float(row["denoised_r1"]))
            per_step["merged"][row["sample"]] += int(float(row["merged"]))
    for row in read_tsv(f"{v2_dir}/05_merged/{rung}/{marker}/track.tsv"):
        per_step["nonchim"][row["sample"]] = int(float(row["nonchim"]))
    return per_step


def mock_summary(features, expected):
    """Mock metrics from features; expected = {sequence: organism}."""
    sample_sums = Counter()
    for _, _, _, mock in features:
        sample_sums.update(mock)
    mocks = {s for s, n in sample_sums.items() if n >= min_mock}
    reads = {(fid, seq): sum(n for s, n in mock.items() if s in mocks) for fid, seq, _, mock in features}
    total = sum(reads.values())
    exact, contained, organism_reads, other = set(), set(), Counter(), []
    for (fid, seq), n in reads.items():
        if n == 0:
            continue
        hit = expected.get(seq)
        if hit:
            exact.add(seq)
            organism_reads[hit] += n
            continue
        inside = [e for e in expected if e in seq]
        if inside:
            contained.add(inside[0])
            organism_reads[expected[inside[0]]] += n
        elif n >= 0.001 * total:
            other.append(n)
    exact_reads = sum(n for (fid, seq), n in reads.items() if seq in expected)
    return {
        "mocks_used": len(mocks),
        "expected_found_exact": f"{len(exact)} of {len(expected)}",
        "expected_found_exact_or_contained": f"{len(exact | contained)} of {len(expected)}",
        "pct_mock_reads_expected_exact": f"{100 * exact_reads / total:.2f}" if total else "-",
        "pct_mock_reads_expected_exact_or_contained": f"{100 * sum(organism_reads.values()) / total:.2f}" if total else "-",
        "other_asvs_ge_0.1pct": len(other),
        "pct_mock_reads_other": f"{100 * (total - sum(organism_reads.values())) / total:.2f}" if total else "-",
        "organisms": {o: 100 * n / total for o, n in organism_reads.items()} if total else {},
    }


def genus_shares(features, taxonomy):
    """% of biosample target reads per genus, and % of target reads with a genus."""
    target, genera = 0, Counter()
    for fid, _, bio, _ in features:
        genus, is_target = taxonomy.get(fid, (None, False))
        if is_target:
            target += bio
            if genus:
                genera[genus] += bio
    return {g: 100 * n / target for g, n in genera.items()}, 100 * sum(genera.values()) / target


remnants = {(row["source"], row["marker"]): row for row in read_tsv(f"{v2_dir}/qc/05_primer_remnants.tsv")}
results = {}
for marker, out_path in [("16S", out_16S), ("ITS1", out_ITS1)]:
    expected = {seq: name.split("|")[0] for name, seq in read_fasta(f"{mock_dir}/{MOCK_FILES[marker]}").items()}
    feats = {"v1": features_from_biom(f"{v1_dir}/export/CFM_{marker}_dada2_table/feature-table.biom",
                                      f"{v1_dir}/export/CFM_{marker}_dada2_repseqs/dna-sequences.fasta")}
    for rung in RUNGS:
        feats[rung] = features_from_seqtab(f"{v2_dir}/05_merged/{rung}/{marker}/seqtab_nochim.tsv")
    steps = {"v1": v1_stats(marker), **{rung: v2_steps(marker, rung) for rung in RUNGS}}

    rows = []   # (section, metric, {dataset: value})
    for scope in ["all", "biosamples"]:
        sums = {d: {step: sum(n for s, n in steps[d][step].items() if scope == "all" or is_biosample(s))
                    for step in STEPS} for d in DATASETS}
        for step in STEPS:
            rows.append(("retention", f"{scope}_{step}_pairs", {d: sums[d][step] for d in DATASETS}))
        for step in STEPS[1:]:
            rows.append(("retention", f"{scope}_{step}_pct_of_demultiplexed",
                         {d: f"{100 * sums[d][step] / sums[d]['demultiplexed']:.2f}" for d in DATASETS}))
    rows.append(("retention", "biosample_median_pct_nonchim", {
        d: f"{median(100 * steps[d]['nonchim'].get(s, 0) / n for s, n in steps[d]['demultiplexed'].items() if is_biosample(s) and n):.2f}"
        for d in DATASETS}))

    rows.append(("features", "asvs", {d: len(feats[d]) for d in DATASETS}))
    rows.append(("features", "asvs_in_biosamples", {d: sum(1 for f in feats[d] if f[2] > 0) for d in DATASETS}))
    rows.append(("features", "median_asv_length", {d: median(len(f[1]) for f in feats[d]) for d in DATASETS}))
    if marker == "ITS1":
        rows.append(("features", "otu97", {
            "v1": len(biom.load_table(f"{v1_dir}/export/CFM_ITS1_otu97_table/feature-table.biom").ids(axis="observation")),
            "C": len(biom.load_table(f"{v2_dir}/export/CFM_ITS1_otu97_table/feature-table.biom").ids(axis="observation"))}))

    for metric, column in [("pct_reads_any_primer_remnant", "pct_reads_any_remnant"), ("asv_length_p01", "len_p01")]:
        rows.append(("remnants", metric, {d: remnants[("v1" if d == "v1" else f"v2_{d}", marker)][column]
                                          for d in DATASETS}))

    mock = {d: mock_summary(feats[d], expected) for d in DATASETS}
    for metric in ["mocks_used", "expected_found_exact", "expected_found_exact_or_contained",
                   "pct_mock_reads_expected_exact", "pct_mock_reads_expected_exact_or_contained",
                   "other_asvs_ge_0.1pct", "pct_mock_reads_other"]:
        rows.append(("mock", metric, {d: mock[d][metric] for d in DATASETS}))
    for organism in sorted(set(expected.values())):
        rows.append(("mock", f"pct_{organism}", {d: f"{mock[d]['organisms'].get(organism, 0):.2f}" for d in DATASETS}))

    # Top genera: v1 / v1_reclassified / C at the level of the main data
    if marker == "16S":
        taxa_feats = {"v1": feats["v1"], "v1_reclassified": feats["v1"], "C": feats["C"]}
        taxonomies = {"v1": f"{v1_dir}/export/CFM_16S_taxonomy_diverse_weighted/taxonomy.tsv",
                      "v1_reclassified": f"{v2_dir}/taxonomy/CFM_16S_v1_taxonomy_{main_16S}.tsv",
                      "C": f"{v2_dir}/taxonomy/CFM_16S_taxonomy_{main_16S}.tsv"}
    else:
        v1_otus = features_from_biom(f"{v1_dir}/export/CFM_ITS1_otu97_table/feature-table.biom")
        taxa_feats = {"v1": v1_otus, "v1_reclassified": v1_otus,
                      "C": features_from_biom(f"{v2_dir}/export/CFM_ITS1_otu97_table/feature-table.biom")}
        taxonomies = {"v1": f"{v1_dir}/export/CFM_ITS1_otu97_taxonomy/taxonomy.tsv",
                      "v1_reclassified": f"{v2_dir}/taxonomy/CFM_ITS1_v1_otu97_taxonomy_{main_ITS1}.tsv",
                      "C": f"{v2_dir}/taxonomy/CFM_ITS1_otu97_taxonomy_{main_ITS1}.tsv"}
    shares, genus_pct = {}, {}
    for d, path in taxonomies.items():
        shares[d], genus_pct[d] = genus_shares(taxa_feats[d], read_taxonomy(path, marker))
    rows.append(("top_genera", "pct_target_reads_with_genus", {d: f"{v:.2f}" for d, v in genus_pct.items()}))
    top10 = {d: [g for g, _ in Counter(s).most_common(10)] for d, s in shares.items()}
    # number of genera two top-10 lists share, in the column of the second dataset
    for a, b in [("v1", "v1_reclassified"), ("v1_reclassified", "C"), ("v1", "C")]:
        rows.append(("top_genera", f"top10_shared_with_{a}", {b: len(set(top10[a]) & set(top10[b]))}))
    union = sorted(set().union(*top10.values()), key=lambda g: -shares["C"].get(g, 0))
    for genus in union:
        rows.append(("top_genera", genus, {d: f"{shares[d].get(genus, 0):.2f}" for d in shares}))

    columns = ["v1", "v1_reclassified"] + list(RUNGS)
    with open(out_path, "w") as out:
        out.write("\t".join(["section", "metric"] + columns) + "\n")
        for section, metric, values in rows:
            out.write("\t".join([section, metric] + [str(values.get(c, "-")) for c in columns]) + "\n")
    results[marker] = {"steps": steps, "mock": mock, "shares": shares, "top10": top10, "union": union}

# ── Figure: rows = markers; columns = retention, mock composition, top genera ──
colours = dict(zip(DATASETS + ["v1_reclassified"], plt.cm.Set2.colors))  # one colour per dataset
fig, axes = plt.subplots(2, 3, figsize=(18, 10))         # overall size in inches
for row, marker in enumerate(["16S", "ITS1"]):
    r = results[marker]
    ax = axes[row, 0]                                       # (1) % of demultiplexed pairs per step
    for d in DATASETS:
        sums = [sum(r["steps"][d][step].values()) for step in STEPS]
        ax.plot(STEPS, [100 * x / sums[0] for x in sums], marker="o", label=d, color=colours[d])
    ax.set_ylabel("% of demultiplexed read pairs")
    ax.set_title(f"{marker}: reads kept per step (all samples)", fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)  # tilted step names
    ax.legend(fontsize=8)

    ax = axes[row, 1]                                       # (2) mock composition, stacked bars
    organisms = sorted({o for d in DATASETS for o in r["mock"][d]["organisms"]})
    org_colours = dict(zip(organisms, plt.cm.tab10.colors))  # one colour per mock organism
    bottom = [0.0] * len(DATASETS)
    for o in organisms + ["other"]:
        heights = [r["mock"][d]["organisms"].get(o, 0) if o != "other"
                   else 100 - sum(r["mock"][d]["organisms"].values()) for d in DATASETS]
        ax.bar(DATASETS, heights, bottom=bottom, label=o.replace("_", " "),
               color=org_colours.get(o, "#bdbdbd"))          # grey = reads outside expected sequences
        bottom = [b + h for b, h in zip(bottom, heights)]
    ax.set_ylabel("% of mock reads")
    ax.set_title(f"{marker}: mock composition", fontsize=10)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3,
              frameon=False)                        # legend below the bars, 3 columns

    ax = axes[row, 2]                                       # (3) top genera, grouped horizontal bars
    genera = r["union"][:15][::-1]                          # at most 15 genera, most abundant in C on top
    height = 0.27                                           # bar thickness per dataset
    for k, d in enumerate(["v1", "v1_reclassified", "C"]):
        ax.barh([i + (k - 1) * height for i in range(len(genera))],
                [r["shares"][d].get(g, 0) for g in genera], height=height, label=d, color=colours[d])
    ax.set_yticks(range(len(genera)))
    ax.set_yticklabels(genera, fontsize=7, style="italic")
    ax.set_xlabel("% of biosample target reads")
    ax.set_title(f"{marker}: top genera (biosamples)", fontsize=10)
    ax.legend(fontsize=8)
fig.suptitle("v1 vs v2 with ablation rungs (A: runs pooled, loess; B: per run, loess; C: per run, binned)",
             fontweight="bold")
fig.tight_layout()
fig.savefig(figure_out, dpi=300)
print(f"Written: {out_16S}, {out_ITS1}, {figure_out}")
