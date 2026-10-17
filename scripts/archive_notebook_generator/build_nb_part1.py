"""Build notebooks/03_16S_ITS1_RFC.ipynb - part 1: cells 0..13 (setup -> helpers)."""
import json, pathlib

CELLS = []
def md(src):   CELLS.append(('markdown', src))
def code(src): CELLS.append(('code', src))

# ── 0. Title ─────────────────────────────────────────────────────────────────
md("""# 03 — Random Forest: 16S ASV & ITS1 OTU97

CLR-transformed amplicon abundances as features. Every target is modelled on three
feature sets in parallel: **ITS1 OTU97**, **16S ASV**, and **Combined** (both concatenated).

| # | Target | Type | Samples | Question |
|---|--------|------|---------|----------|
| 1 | `farm_id` | 7 classes | all flowers | does the microbiome carry a site signature? |
| 2 | `sample_type` | binary | all flowers | bagged vs unbagged — visiting-insect footprint |
| 3 | `pollination` | 3 classes | unbagged only | none / weak / strong — pollinator footprint |
| 4 | `fungal_colonization` | binary | unbagged only | stigma colonisation score 0 vs 5 |
| 5 | `pi_clsm` | **regression** | unbagged only | pollen load as a continuous response |

`management_type` is deliberately **not** a target: each farm has exactly one management
type (7 farms x 42 flowers, no overlap), so a model separating management levels is
separating farms. It would duplicate target 1 under a different name.

**Evaluation.** Hyperparameters are tuned by nested cross-validation: an inner
`GridSearchCV` selects them inside every outer fold, so the reported score is not
inflated by the selection itself. Statistical significance comes from a permutation
test against shuffled labels, not from a fixed chance line.

Metrics: balanced accuracy, MCC, ROC-AUC (OvR macro), permutation p-value, confusion
matrix, calibration; Gini / permutation / SHAP importance.""")

# ── 1. Imports & parameters ──────────────────────────────────────────────────
md("## 0. Imports & parameters")

code('''import warnings; warnings.filterwarnings('ignore')
import json, os, time
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import biom
from scipy.stats import spearmanr

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import (StratifiedKFold, RepeatedStratifiedKFold,
                                     KFold, RepeatedKFold, GridSearchCV,
                                     cross_validate, cross_val_predict,
                                     permutation_test_score)
from sklearn.inspection import permutation_importance as sklearn_perm_imp
from sklearn.metrics import (confusion_matrix, ConfusionMatrixDisplay,
                             balanced_accuracy_score, matthews_corrcoef,
                             roc_auc_score, brier_score_loss, make_scorer,
                             r2_score, mean_absolute_error)
from sklearn.calibration import calibration_curve
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.pipeline import Pipeline
import shap

ROOT      = Path('..').resolve()
BIOM_DIR  = ROOT / 'results' / 'biom'
FIG_DIR   = ROOT / 'results' / 'figures' / 'random_forest'
TAB_DIR   = ROOT / 'results' / 'tables'
CACHE_DIR = TAB_DIR / 'rf_cache'
for d in (FIG_DIR, TAB_DIR, CACHE_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ── Output conventions ───────────────────────────────────────────────────────
# RFC_TAG marks where a result was produced, so outputs from the cluster and from
# the laptop can sit in the same directory without overwriting each other after an
# rsync. Set RFC_TAG=_hpc in the SLURM script; empty locally.
OUT_TAG    = os.environ.get('RFC_TAG', '')
FIG_PREFIX = f'03{OUT_TAG}_'   # matches the notebook number; old 06a_* outputs are kept
FIG_DPI    = 300               # project rule: every figure is saved at 300 dpi

# ── Run profile ──────────────────────────────────────────────────────────────
# QUICK=True  : fast pass for developing / debugging the notebook (~25-30 min)
# QUICK=False : full run whose numbers go into the manuscript (~4 h on 6 cores)
QUICK           = True
FORCE_RECOMPUTE = False   # True = ignore the cache and recompute everything

RANDOM_STATE = 42
N_SPLITS     = 5          # folds in every CV loop (outer and OOF)
N_TOP_IMP    = 20         # bars shown in importance plots

if QUICK:
    N_TREES       = 200
    N_REPEATS     = 1     # repeats of the outer CV loop
    INNER_SPLITS  = 3
    N_PERM        = 50    # label shuffles, ordinary targets
    N_PERM_STRONG = 100   # label shuffles, targets whose p sits on the resolution floor
    PARAM_GRID    = {'min_samples_leaf': [1, 5], 'max_depth': [None, 10]}
else:
    N_TREES       = 500
    N_REPEATS     = 3
    INNER_SPLITS  = 5
    N_PERM        = 200
    N_PERM_STRONG = 1000
    # min_samples_leaf=1 and max_depth=None are the sklearn defaults, kept in the
    # grid so tuning is free to select the untuned configuration.
    PARAM_GRID    = {'min_samples_leaf': [1, 2, 5, 10], 'max_depth': [None, 10, 20]}

# Permutation tests use fewer trees than the main models. This is sound: the observed
# score and the null distribution go through the *same* procedure, which is what makes
# the p-value valid. The observed score therefore differs slightly from the nested-CV
# estimate -- a different procedure, not a discrepancy.
N_TREES_PERM = 200

# Repeats per feature inside permutation importance (a different thing from the
# permutation significance test above).
N_PERM_IMP = 5 if QUICK else 10

# Combinations whose p-value hit the 1/(N_PERM+1) floor in the first full run; they
# need more shuffles to clear the Bonferroni threshold (0.05 / 12 tests = 0.0042).
STRONG_COMBOS = {
    ('farm_id', 'ITS1'), ('farm_id', '16S'), ('farm_id', 'Combined'),
    ('sample_type', 'ITS1'), ('sample_type', '16S'), ('sample_type', 'Combined'),
    ('fungal_colonization', 'Combined'),
}

# Trees are built with n_jobs=-1; the CV loops around them stay sequential so the
# two levels of parallelism do not oversubscribe the 6 available cores.
N_JOBS_CV = 1

# ── Feature filtering ────────────────────────────────────────────────────────
# 10% prevalence matches 08a/08b and the ANCOMBC2 scripts (PREV_THRESH <- 0.10).
# Keeping the same feature space is what makes RF importances comparable with
# ANCOMBC2 / ISA results. No abundance filter is applied anywhere in this project.
PREV_THRESH = 0.10

# ── Sample exclusions ────────────────────────────────────────────────────────
# cyb43 : bagged flower with pi_clsm=340 -> the bag failed, the 'bagged' label is wrong.
#         Dropped everywhere.
# cvi11, cvi51 : bagged flowers with an insect found inside during specimen prep.
#         They violate the 'no insect access' definition, so they are dropped from
#         sample_type. They are bagged, so unbagged-only targets exclude them anyway.
#         farm_id keeps them: farm identity is unaffected by insect contact.
OUTLIER_IDS = ['cyb43']
INSECT_IDS  = ['cvi11', 'cvi51']

PI_THRESH_POLL = 40   # pollen grains per stigma; none / weak / strong boundary (as in 07a/07b)

# ── Plot colours ─────────────────────────────────────────────────────────────
# One colour per feature set, reused in every panel so the reader can track them.
DS_COLOURS = {'ITS1': '#9467BD', '16S': '#D62728', 'Combined': '#2CA02C'}

print(f'Profile: {"QUICK" if QUICK else "FULL"}  |  trees={N_TREES}  '
      f'outer={N_SPLITS}x{N_REPEATS}  inner={INNER_SPLITS}  '
      f'grid={len(PARAM_GRID["min_samples_leaf"]) * len(PARAM_GRID["max_depth"])} combos  '
      f'perm={N_PERM}/{N_PERM_STRONG} @ {N_TREES_PERM} trees')
print(f'Figures -> {FIG_DIR}')
print(f'Cache   -> {CACHE_DIR}')''')

# ── 2. Load data ─────────────────────────────────────────────────────────────
md("## 1. Load data")

md("""### 1a. Taxonomy parsing

Rank prefixes (`g__`, `f__`, ...) are **kept**, because the R side of the project keeps
them too: the `tax_table` columns feeding the ANCOMBC2 and ISA figures hold values like
`f__Bulleribasidiaceae` and `g__Hannaella`. Identical labels let an RF importance bar
chart sit next to an ANCOMBC2 volcano plot and be read as one figure.

`deepest_rank` walks down the hierarchy until it finds an informative rank, mirroring
`make_tax_label` in `scripts/16a_16S_ancombc2_insect_footprint.R`. It replaces the old
`unclassified` bin, which pooled 55 of 303 ITS1 features and 47 of 231 16S features into
a single row and consequently topped the genus heatmap.""")

code('''# Values that carry no information even when the rank field is populated.
TAX_SKIP = {'', 'NA', 'nan', 'unidentified', 'uncultured', 'uncultured_bacterium',
            'metagenome', 'Incertae_sedis', 'unidentified_Fungi', 'Fungi_sp'}

RANK_ORDER = ['g', 'f', 'o', 'c', 'p', 'k']   # deepest first, used by the fallback


def _rank_map(taxon_str):
    """Split a taxonomy string into {rank letter: 'g__Name'}, dropping empty ranks.

    UNITE species-hypothesis fields (sh__) are skipped; QIIME leaves a bare 'g__' when
    the rank is unassigned, which is filtered out by the TAX_SKIP check.
    """
    out = {}
    for part in str(taxon_str).split(';'):
        part = part.strip()
        if '__' not in part or part.startswith('sh__'):
            continue
        rank, val = part.split('__', 1)
        val = (val.strip().replace('_gen_Incertae_sedis', '')
                          .replace('_fam_Incertae_sedis', ''))
        if val and val not in TAX_SKIP:
            out[rank.strip().lower()] = f'{rank.strip().lower()}__{val}'
    return out


def deepest_rank(taxon_str, feature_id=''):
    """Deepest informative rank, prefix kept: g__Fusarium -> f__Nectriaceae -> ...

    Used as the aggregation key in the genus-level heatmap.
    """
    ranks = _rank_map(taxon_str)
    for r in RANK_ORDER:
        if r in ranks:
            return ranks[r]
    return str(feature_id)[:12] or 'unassigned'


def feature_label(taxon_str, n_levels=2):
    """Deepest n informative ranks with prefixes, for per-feature plot labels.

    e.g. 'f__Bulleribasidiaceae g__Vishniacozyma'. Longer than the bare names used
    previously, but matches how taxa are labelled in the R figures.
    """
    ranks = _rank_map(taxon_str)
    ordered = [ranks[r] for r in reversed(RANK_ORDER) if r in ranks]   # k -> g
    return ' '.join(ordered[-n_levels:]) if ordered else str(taxon_str)[:20]


def load_marker(biom_stem, label):
    """Load a feature table plus its taxonomy; returns counts and two label maps."""
    counts = biom.load_table(str(BIOM_DIR / f'{biom_stem}_feature-table.biom')) \\
               .to_dataframe(dense=True)
    tax_raw = pd.read_csv(BIOM_DIR / f'{biom_stem}_taxonomy.tsv', sep='\\t', index_col=0)
    labels  = {i: feature_label(r['Taxon'])      for i, r in tax_raw.iterrows()}
    deepest = {i: deepest_rank(r['Taxon'], i)    for i, r in tax_raw.iterrows()}
    print(f'{label}: {counts.shape[0]} features  {counts.shape[1]} samples')
    return counts, labels, deepest


counts_its1, tax_its1, rank_its1 = load_marker(
    'cfm_its1_otu97_fungi_biosamples', 'ITS1 OTU97')
counts_16s, tax_16s, rank_16s = load_marker(
    'cfm_16s_bacteria_biosamples', '16S ASV')

_ex = list(tax_its1.items())[:2]
print('\\nExample labels:', *[f'\\n  {v}' for _, v in _ex])''')

md("### 1c. Metadata")

code('''meta_raw = pd.read_csv(ROOT / 'data' / 'sample_metadata_with_microscopy.txt',
                       sep='\\t', index_col=0, dtype=str)
meta_raw = meta_raw.drop(index='#q2:types', errors='ignore').copy()   # QIIME2 type row
meta_raw.index.name = 'sample_id'

# Decimal commas appear in the microscopy columns; coerce to float.
for col in ['pi_clsm', 'fungal_colonization_score']:
    meta_raw[col] = pd.to_numeric(
        meta_raw[col].astype(str).str.replace(',', '.', regex=False), errors='coerce')

meta_raw = meta_raw[meta_raw['sample_type'].isin(['bagged_flower', 'unbagged_flower'])]
n_all = len(meta_raw)

meta = meta_raw[~meta_raw.index.isin(OUTLIER_IDS)].copy()

print(f'{n_all:4d}  all bagged + unbagged flowers')
print(f'{len(meta):4d}  after dropping {OUTLIER_IDS} (bag failure)')
print(f'{len(meta[~meta.index.isin(INSECT_IDS)]):4d}  after also dropping {INSECT_IDS} '
      f'(insect inside bag) -> used for sample_type')
print()
print(meta['sample_type'].value_counts().to_string())''')

# ── 3. Transformations & filtering ───────────────────────────────────────────
md("""## 2. Feature filtering and transformation

**Prevalence filter (10%) and the transform are computed on all 294 flowers, before
any target-specific subsetting.** That is safe for CLR: CLR is a per-sample (per-row)
operation, so selecting a subset of samples afterwards does not change the feature
values. It also keeps every target in the same feature space, and the same space that
08a/08b and the ANCOMBC2 scripts use.

Four transforms are implemented. CLR is the default; the other three are used by the
sensitivity panel further down.

| Transform | What it does | Note |
|-----------|--------------|------|
| `raw` | untransformed counts | library size stays in the features |
| `tss` | counts / library size | removes depth, keeps the compositional constraint |
| `clr` | multiplicative zero replacement, then centred log-ratio | **default** |
| `rclr` | log of non-zero values minus the row mean of those logs; zeros stay 0 | matches `vegan::vegdist(method='robust.aitchison')` used for the PCoA in 08a/09c |

Random Forest is *not* invariant to these choices. Trees are invariant to monotonic
transforms applied per feature, but CLR and TSS divide by a per-sample quantity, which
can reorder samples within a feature — and that ordering is exactly what a split uses.""")

code('''def multiplicative_replacement(mat, delta=None):
    """Replace zeros with a small constant, rescaling the rest to keep row sums.

    mat: samples x features, rows are relative abundances.
    Simple fixed-delta variant (skbio default). The R pipeline (08a) uses
    zCompositions::cmultRepl(method='CZM'), which estimates the replacement per
    feature; values therefore differ slightly between the two implementations.
    """
    mat = np.asarray(mat, dtype=float)
    if delta is None:
        delta = 1.0 / mat.shape[1] ** 2
    out = mat.copy()
    for i in range(len(out)):
        zm = out[i] == 0
        k  = zm.sum()
        if k:
            out[i,  zm] = delta
            out[i, ~zm] *= 1.0 - k * delta
    return out


def transform_counts(counts_f, how='clr'):
    """counts_f: features x samples -> DataFrame samples x features."""
    ra = counts_f.div(counts_f.sum(axis=0), axis=1)      # relative abundance
    idx, cols = counts_f.columns, counts_f.index

    if how == 'raw':
        return pd.DataFrame(counts_f.T.values, index=idx, columns=cols)
    if how == 'tss':
        return pd.DataFrame(ra.T.values, index=idx, columns=cols)
    if how == 'clr':
        comp = multiplicative_replacement(ra.T.values)
        logc = np.log(comp)
        return pd.DataFrame(logc - logc.mean(axis=1, keepdims=True), index=idx, columns=cols)
    if how == 'rclr':
        # Robust CLR, as implemented by vegan: the geometric mean is taken over the
        # non-zero entries of each sample and zeros are left at 0.
        X    = ra.T.values
        nz   = X > 0
        logX = np.where(nz, np.log(np.where(nz, X, 1.0)), 0.0)
        row_gm = logX.sum(axis=1) / np.maximum(nz.sum(axis=1), 1)
        return pd.DataFrame(np.where(nz, logX - row_gm[:, None], 0.0), index=idx, columns=cols)
    raise ValueError(f'unknown transform: {how}')


def prevalence_filter(counts_raw, thresh, subset_ids=None):
    """Keep features present in >= thresh of samples.

    subset_ids: compute prevalence within this subset instead of all samples
    (used by the sensitivity panel to test the unbagged-only variant).
    """
    ref = counts_raw if subset_ids is None else counts_raw[
        [s for s in counts_raw.columns if s in set(subset_ids)]]
    prev = (ref > 0).sum(axis=1) / ref.shape[1]
    return counts_raw[prev >= thresh]


def build_features(counts_raw, meta_df, tax_dict, how='clr', thresh=PREV_THRESH,
                   subset_ids=None, label=''):
    """Filter -> transform -> align to metadata. Returns (X, feature labels)."""
    counts_f = prevalence_filter(counts_raw, thresh, subset_ids)
    X = transform_counts(counts_f, how)
    X = X.loc[X.index.intersection(meta_df.index)]
    feats = [f'{tax_dict.get(f, str(f)[:18])} [{str(f)[:6]}]' for f in X.columns]
    if label:
        print(f'{label:10s} {how:5s}: {X.shape[1]:5d} features  {X.shape[0]:3d} samples')
    return X, feats''')

code('''X_its1, feat_its1 = build_features(counts_its1, meta, tax_its1, 'clr', label='ITS1 OTU97')
X_16s,  feat_16s  = build_features(counts_16s,  meta, tax_16s,  'clr', label='16S ASV')

# Combined multi-kingdom matrix: samples present in both markers, columns prefixed
# so ITS1 and 16S feature names cannot collide.
common_both = X_its1.index.intersection(X_16s.index)
X_comb = pd.concat(
    [X_its1.loc[common_both].add_prefix('ITS1_'), X_16s.loc[common_both].add_prefix('S16_')],
    axis=1)
feat_comb = ([f'ITS1 {f}' for f in feat_its1] + [f'16S {f}' for f in feat_16s])
print(f'{"Combined":10s} clr  : {X_comb.shape[1]:5d} features  {X_comb.shape[0]:3d} samples')

# Farm code per sample, aligned to each feature matrix. Needed by the within-farm
# centering transformer, which receives it as an extra column of X.
farm_codes = pd.Categorical(meta['farm_id']).codes
FARM = pd.Series(farm_codes, index=meta.index, name='farm')

# Registry consumed by the target x dataset loop below.
# 'rank' maps feature id -> deepest informative rank (g__Fusarium, f__Nectriaceae, ...)
# and is the aggregation key of the taxonomic importance heatmap.
FEATURES = {
    'ITS1':     dict(X=X_its1, feats=feat_its1, ids=list(X_its1.columns), rank=rank_its1),
    '16S':      dict(X=X_16s,  feats=feat_16s,  ids=list(X_16s.columns),  rank=rank_16s),
    'Combined': dict(X=X_comb, feats=feat_comb, ids=list(X_comb.columns), rank=None),
}
DATASETS = list(FEATURES)''')

pathlib.Path(__file__).with_name('cells_part1.json').write_text(
    json.dumps(CELLS, ensure_ascii=False))
print(f'part 1: {len(CELLS)} cells')
