#!/usr/bin/env python3
"""
Heatmap of log10 relative abundances for ITS1 fungal ASVs (top-N by mean RA).

Inputs (defaults match current repo):
  - BIOM table: results/biom/cfm_its1_fungi_biosamples_feature-table.biom
  - Taxonomy TSV: results/biom/cfm_its1_fungi_biosamples_taxonomy.tsv
  - QIIME2 metadata: data/cfm_qiime2_metadata.txt

Output:
  - PNG heatmap in results/figures/
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from biom import load_table
from matplotlib import pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--biom",
        default="results/biom/cfm_its1_fungi_biosamples_feature-table.biom",
        help="BIOM feature table (ASV x sample).",
    )
    p.add_argument(
        "--taxonomy",
        default="results/biom/cfm_its1_fungi_biosamples_taxonomy.tsv",
        help="TSV with columns: Feature ID, Taxon.",
    )
    p.add_argument(
        "--metadata",
        default="data/cfm_qiime2_metadata.txt",
        help="QIIME2 metadata TSV (2nd row '#q2:types' is ignored).",
    )
    p.add_argument("--top-n", type=int, default=50, help="Top N ASVs by mean relative abundance.")
    p.add_argument("--eps", type=float, default=1e-6, help="Pseudocount for log10(RA + eps).")
    p.add_argument(
        "--out",
        default="results/figures/ITS1_ASV_log10RA_heatmap_top{top_n}.png",
        help="Output PNG path (can include {top_n}).",
    )
    return p.parse_args()


def taxon_to_label(taxon: str) -> str:
    # Prefer species, then genus, then last available rank.
    parts = [p.strip() for p in taxon.split(";") if p.strip()]
    if not parts:
        return "unclassified"

    def find_prefix(prefix: str) -> str | None:
        for p in reversed(parts):
            if p.startswith(prefix) and p != prefix:
                return p
        return None

    species = find_prefix("s__")
    genus = find_prefix("g__")
    chosen = species or genus or parts[-1]
    chosen = re.sub(r"^[a-z]__", "", chosen)
    chosen = chosen.replace("_", " ")
    return chosen or "unclassified"


def make_unique(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
        else:
            seen[n] += 1
            out.append(f"{n}_{seen[n]}")
    return out


def main() -> None:
    args = parse_args()

    biom_fp = Path(args.biom)
    tax_fp = Path(args.taxonomy)
    meta_fp = Path(args.metadata)

    if not biom_fp.exists():
        raise SystemExit(f"Missing BIOM: {biom_fp}")
    if not tax_fp.exists():
        raise SystemExit(f"Missing taxonomy TSV: {tax_fp}")
    if not meta_fp.exists():
        raise SystemExit(f"Missing metadata: {meta_fp}")

    table = load_table(str(biom_fp))
    # observations (ASVs) x samples
    df_counts: pd.DataFrame = table.to_dataframe(dense=True)

    # Relative abundance per sample (TSS)
    col_sums = df_counts.sum(axis=0)
    df_ra = df_counts.div(col_sums.replace(0, np.nan), axis=1).fillna(0.0)

    top_n = int(args.top_n)
    top_asvs = df_ra.mean(axis=1).sort_values(ascending=False).head(top_n).index.tolist()
    df_top = df_ra.loc[top_asvs]

    df_log = np.log10(df_top + float(args.eps))

    # Taxonomy labels
    tax = pd.read_csv(tax_fp, sep="\t")
    feature_col = "Feature ID" if "Feature ID" in tax.columns else tax.columns[0]
    tax = tax.set_index(feature_col)
    if "Taxon" in tax.columns:
        labels = tax["Taxon"].reindex(df_log.index).fillna("unclassified").map(taxon_to_label)
    else:
        labels = pd.Series(index=df_log.index, data=df_log.index)
    # ensure uniqueness
    df_log.index = make_unique(labels.astype(str).tolist())

    # Metadata for column color bars (optional)
    meta = pd.read_csv(meta_fp, sep="\t", comment="#")
    if "sample-id" in meta.columns:
        meta = meta.set_index("sample-id")
    else:
        meta = meta.set_index(meta.columns[0])

    meta = meta.reindex(df_log.columns)
    col_colors = []

    def add_col_color(col: str, palette: str = "tab10") -> None:
        if col not in meta.columns:
            return
        vals = meta[col].astype("string")
        cats = [c for c in sorted(vals.dropna().unique().tolist()) if c != "<NA>"]
        pal = sns.color_palette(palette, n_colors=max(3, len(cats))).as_hex()
        lut = dict(zip(cats, pal[: len(cats)]))
        col_colors.append(vals.map(lut).fillna("#D9D9D9"))

    add_col_color("sample_type", "Set2")
    add_col_color("management_type", "Set3")
    add_col_color("farm_id", "tab20")

    col_colors_df = None
    if col_colors:
        col_colors_df = pd.concat(col_colors, axis=1)
        col_colors_df.columns = [c for c in ["sample_type", "management_type", "farm_id"] if c in meta.columns]

    out_fp = Path(args.out.format(top_n=top_n))
    out_fp.parent.mkdir(parents=True, exist_ok=True)

    sns.set_context("notebook")
    sns.set_style("white")

    # Figure size heuristic: scale with number of samples a bit, but cap.
    n_samples = df_log.shape[1]
    fig_w = min(30, max(10, n_samples / 18))
    fig_h = 10

    g = sns.clustermap(
        df_log,
        method="average",
        metric="euclidean",
        cmap="mako",
        col_colors=col_colors_df,
        xticklabels=False,
        yticklabels=True,
        figsize=(fig_w, fig_h),
        cbar_kws={"label": f"log10(RA + {args.eps:g})"},
    )

    title = f"ITS1 fungi: top {top_n} ASVs by mean relative abundance\nValues: log10(RA + eps)"
    g.fig.suptitle(title, y=1.02, fontsize=12)
    g.ax_heatmap.set_xlabel("Samples")
    g.ax_heatmap.set_ylabel("ASVs (taxonomy label)")

    g.fig.savefig(out_fp, dpi=200, bbox_inches="tight")
    plt.close(g.fig)

    print(f"Wrote: {out_fp}")


if __name__ == "__main__":
    main()
