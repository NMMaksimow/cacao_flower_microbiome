#!/usr/bin/env python3
"""Random Forest analysis of the cacao flower microbiome -- batch version.

Generated from notebooks/03_16S_ITS1_RFC.ipynb; edit the notebook, not this file.

Run from the notebooks/ directory, which is what the relative paths assume:
    cd notebooks && python 03_16S_ITS1_RFC.py

Environment variables:
    RFC_QUICK=1   fast pass for checking the environment (minutes, numbers not usable)
    RFC_QUICK=0   full run for the manuscript (default)
"""
import os
import matplotlib
matplotlib.use('Agg')          # no display on a compute node; must precede pyplot


# ===== cell 0 ==========================================================
import warnings; warnings.filterwarnings('ignore')
import hashlib, json, os, re, time
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
# RFC_TAG marks where a result was produced, so cluster and laptop outputs can share
# a directory after an rsync. It goes after the script number, keeping files sorted by
# step: RFC_TAG=_hpc -> 03_hpc_summary_heatmap.png ; unset -> 03_summary_heatmap.png
OUT_TAG    = os.environ.get('RFC_TAG', '')
FIG_PREFIX = f'03{OUT_TAG}_'   # matches the notebook number; old outputs are kept
FIG_DPI    = 300               # project rule: every figure is saved at 300 dpi

# ── Run profile ──────────────────────────────────────────────────────────────
# QUICK=True  : fast pass for developing / debugging the notebook (~25-30 min)
# QUICK=False : full run whose numbers go into the manuscript (~4 h on 6 cores)
QUICK = os.environ.get('RFC_QUICK', '0') == '1'
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
# (target, dataset, variant) triples. The within-farm sample_type rows were added
# after the first full run: they carry the paper's main result, and at 200 shuffles
# their p sat on the 1/(n+1) floor of 0.005, which does not clear the Bonferroni
# threshold of 0.05/27 = 0.00185. At 1000 the floor drops to 0.001 and it does.
STRONG_COMBOS = {
    ('farm_id', 'ITS1', 'raw'), ('farm_id', '16S', 'raw'),
    ('farm_id', 'Combined', 'raw'),
    ('sample_type', 'ITS1', 'raw'), ('sample_type', '16S', 'raw'),
    ('sample_type', 'Combined', 'raw'),
    ('fungal_colonization', 'Combined', 'raw'),
    ('sample_type', 'ITS1', 'within-farm'), ('sample_type', '16S', 'within-farm'),
    ('sample_type', 'Combined', 'within-farm'),
}

# Parallelism goes on the CV loops, not on tree building.
#
# Measured on a 32-core node: with n_jobs=-1 on the forest and sequential CV loops, a
# run reported real 109m10s against user 107m56s -- a ratio of 0.99, i.e. one core busy
# out of 32. The samples are small (293 x 303), so splitting 500 trees across workers
# costs more in dispatch than it saves. The outer loops are where the work is: 915
# independent fits in nested CV, n_permutations x 5 folds in the significance test,
# n_features x n_repeats x 5 folds in permutation importance.
#
# Results are unaffected: a forest with a fixed random_state is deterministic
# regardless of thread count, and the CV splits are seeded too.
N_JOBS_CV = -1

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
print(f'Cache   -> {CACHE_DIR}')

# ===== cell 1 ==========================================================
# Values that carry no information even when the rank field is populated.
TAX_SKIP = {'', 'NA', 'nan', 'unidentified', 'uncultured', 'uncultured_bacterium',
            'metagenome', 'Incertae_sedis', 'Incertae_Sedis', 'Unknown_Family',
            'unidentified_Fungi', 'Fungi_sp'}

# Two top ranks, because the two references disagree on the word: SILVA labels the
# top rank d__Bacteria, UNITE labels it k__Fungi. No string carries both, so their
# order relative to each other never matters -- 912 of 25035 16S ASVs are classified
# no further than the domain, and without 'd' they would get no rank at all.
RANK_ORDER = ['g', 'f', 'o', 'c', 'p', 'k', 'd']   # deepest first, used by the fallback
RANK_DOWN  = ['d', 'k', 'p', 'c', 'o', 'f', 'g']   # same ranks, top rank first

# An unplaced rank is named after the taxon enclosing it. UNITE writes
# '<taxon>_<rank>_Incertae_sedis' (g__Fungi_gen_Incertae_sedis), SILVA writes
# '<taxon>_Incertae_Sedis' with no rank token and a capital S
# (f__Rhizobiales_Incertae_Sedis). Both reduce to the enclosing taxon, which the
# de-duplication below then recognises as a repeat of the rank above.
_INCERTAE = re.compile(r'_(?:(?:gen|fam|ord|cls|phy)_)?Incertae_[Ss]edis$')


def _rank_map(taxon_str):
    """Split a taxonomy string into {rank letter: 'g__Name'}, dropping ranks that
    carry no information of their own.

    The subtle case is UNITE's unplaced ranks. A string like

        k__Fungi; p__Fungi_phy_Incertae_sedis; ... ; g__Fungi_gen_Incertae_sedis

    means "a fungus, nothing below kingdom is known". Stripping the suffix alone
    invents a genus called Fungi -- 177 OTUs acquired one that way, and the same
    mechanism produced g__Ascomycota (61), g__Agaricomycetes (36), g__Xylariales
    (15) and g__Polyporaceae (7), together 27% of the importance mass in the ITS1
    heatmap. So ranks are walked top-down and a value repeating a higher rank is
    dropped: it marks placement inside that taxon, not a name at this level.

    SILVA does the same thing with f__Rhizobiales_Incertae_Sedis and
    o__Gammaproteobacteria_Incertae_Sedis, and additionally uses the literal
    placeholders 'uncultured' (3331 genera, 1100 families) and 'Unknown_Family';
    those are handled by TAX_SKIP.

    sh__ (species hypothesis) is skipped; QIIME's bare 'g__' is caught by TAX_SKIP.
    """
    raw = {}
    for part in str(taxon_str).split(';'):
        part = part.strip()
        if '__' not in part or part.startswith('sh__'):
            continue
        rank, val = part.split('__', 1)
        raw[rank.strip().lower()] = _INCERTAE.sub('', val.strip())

    out, seen = {}, set()
    for r in RANK_DOWN:                    # top-down, so 'seen' holds higher ranks only
        val = raw.get(r, '')
        if val and val not in TAX_SKIP and val not in seen:
            out[r] = f'{r}__{val}'
            seen.add(val)
    return out


def deepest_rank(taxon_str, feature_id=''):
    """Deepest informative rank, prefix kept: g__Fusarium -> f__Nectriaceae -> ...

    The single label used everywhere: as the aggregation key of the rank heatmap and,
    with the feature hash appended by build_features, as the per-feature plot label.
    One rank rather than two matches make_tax_label() in the ANCOMBC2 / ISA scripts,
    so an RF importance bar chart can sit next to an ANCOMBC2 volcano plot.

    s__ and sh__ are deliberately not candidates: species calls from a short amplicon
    are unreliable, and the R side stops at genus for the same reason.
    """
    ranks = _rank_map(taxon_str)
    for r in RANK_ORDER:
        if r in ranks:
            return ranks[r]
    return str(feature_id)[:12] or 'unassigned'


def load_marker(biom_stem, label):
    """Load a feature table plus its taxonomy; returns counts and one label map."""
    counts = biom.load_table(str(BIOM_DIR / f'{biom_stem}_feature-table.biom')) \
               .to_dataframe(dense=True)
    tax_raw = pd.read_csv(BIOM_DIR / f'{biom_stem}_taxonomy.tsv', sep='\t', index_col=0)
    tax = {i: deepest_rank(r['Taxon'], i) for i, r in tax_raw.iterrows()}
    print(f'{label}: {counts.shape[0]} features  {counts.shape[1]} samples')
    return counts, tax


counts_its1, tax_its1 = load_marker('cfm_its1_otu97_fungi_biosamples', 'ITS1 OTU97')
counts_16s, tax_16s = load_marker('cfm_16s_bacteria_biosamples', '16S ASV')

_ex = list(tax_its1.items())[:2]
print('\nExample labels:', *[f'\n  {v}' for _, v in _ex])

# ===== cell 2 ==========================================================
meta_raw = pd.read_csv(ROOT / 'data' / 'sample_metadata_with_microscopy.txt',
                       sep='\t', index_col=0, dtype=str)
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
print(meta['sample_type'].value_counts().to_string())

# ===== cell 3 ==========================================================
def multiplicative_replacement(mat, delta=None):
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
    return X, feats

# ===== cell 4 ==========================================================
X_its1, feat_its1 = build_features(counts_its1, meta, tax_its1, 'clr', label='ITS1 OTU97')
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
# 'rank' maps feature id -> deepest informative rank (g__Fusarium, f__Nectriaceae, ...),
# the aggregation key of the taxonomic importance heatmap. 'feats' is the same value
# with the feature hash appended, which is what the per-feature plots label bars with.
FEATURES = {
    'ITS1':     dict(X=X_its1, feats=feat_its1, ids=list(X_its1.columns), rank=tax_its1),
    '16S':      dict(X=X_16s,  feats=feat_16s,  ids=list(X_16s.columns),  rank=tax_16s),
    'Combined': dict(X=X_comb, feats=feat_comb, ids=list(X_comb.columns), rank=None),
}
DATASETS = list(FEATURES)

# id -> plot label, applied when a figure is drawn rather than when importances are
# computed, so a cached importance table never carries a stale taxonomy string. Kept
# separate from 'rank' because Combined has no rank map of its own -- its columns
# carry marker prefixes.
for _d in FEATURES.values():
    _d['label'] = dict(zip(_d['ids'], _d['feats']))

# ===== cell 5 ==========================================================
def target_farm_id(m):
    """7 farms, all flowers. Insect-carrying flowers are kept: farm identity is
    unaffected by whether an insect reached the flower."""
    return m['farm_id'].dropna()


def target_sample_type(m):
    """Bagged vs unbagged. cvi11/cvi51 are excluded here: an insect was found inside
    those bags, so they break the class definition this target is testing."""
    return m.loc[~m.index.isin(INSECT_IDS), 'sample_type'].dropna()


def _unbagged(m):
    return m[m['sample_type'] == 'unbagged_flower']


def target_pollination(m):
    """3 classes by pollen load. Thresholds match 07a/07b."""
    s = _unbagged(m)['pi_clsm'].dropna()
    return s.apply(lambda v: 'none' if v == 0 else ('weak' if v < PI_THRESH_POLL else 'strong'))


def target_fungal_colonization(m):
    """Extreme classes of the stigma colonisation score only (0 vs 5).
    Intermediate scores 1-4 are dropped: the contrast is cleaner and classes 3 and 4
    hold 7-8 samples each, too few to cross-validate."""
    s = _unbagged(m)['fungal_colonization_score'].dropna()
    s = s[s.isin([0, 5])]
    return s.map({0.0: 'none', 5.0: 'strong'})


def target_pi_clsm(m):
    """Continuous pollen load, unbagged flowers. Regression target."""
    return _unbagged(m)['pi_clsm'].dropna()


TARGETS = {
    'farm_id':             dict(fn=target_farm_id,             kind='clf', order=None),
    'sample_type':         dict(fn=target_sample_type,         kind='clf',
                                order=['bagged_flower', 'unbagged_flower']),
    'pollination':         dict(fn=target_pollination,         kind='clf',
                                order=['none', 'weak', 'strong']),
    'fungal_colonization': dict(fn=target_fungal_colonization, kind='clf',
                                order=['none', 'strong']),
    'pi_clsm':             dict(fn=target_pi_clsm,             kind='reg', order=None),
}
CLF_TARGETS = [t for t, c in TARGETS.items() if c['kind'] == 'clf']


def build_xy(target, dataset, X_override=None):
    """Align a target vector with a feature matrix. X_override lets the sensitivity
    panels reuse this with a differently transformed / filtered matrix."""
    X = FEATURES[dataset]['X'] if X_override is None else X_override
    y = TARGETS[target]['fn'](meta)
    common = X.index.intersection(y.index)
    return X.loc[common], y.loc[common]


# Sample counts per target, for the record.
rows = []
for t in TARGETS:
    for ds in DATASETS:
        _X, _y = build_xy(t, ds)
        rows.append(dict(target=t, dataset=ds, n=len(_y),
                         classes='—' if TARGETS[t]['kind'] == 'reg'
                                 else str(dict(sorted(_y.value_counts().items())))))
counts_df = pd.DataFrame(rows)
print(counts_df.to_string(index=False))

# ===== cell 6 ==========================================================
def make_clf(**overrides):
    # n_jobs=1: the forest stays single-threaded; parallelism lives in the CV loops
    # around it (see N_JOBS_CV). Nesting both oversubscribes the machine.
    params = dict(n_estimators=N_TREES, max_features='sqrt', class_weight='balanced',
                  random_state=RANDOM_STATE, n_jobs=1)
    params.update(overrides)
    return RandomForestClassifier(**params)


def make_reg(**overrides):
    params = dict(n_estimators=N_TREES, max_features='sqrt',
                  random_state=RANDOM_STATE, n_jobs=1)
    params.update(overrides)
    return RandomForestRegressor(**params)


# ── Within-farm centering ────────────────────────────────────────────────────
class WithinFarmCentering(BaseEstimator, TransformerMixin):
    """Subtract each farm's mean from every feature, removing between-farm variation.

    The farm code arrives as the LAST column of X and is stripped on output. That
    detour exists so the transformer can live inside a CV pipeline: farm means are
    then learned on the training fold only. Centering once outside the CV would let
    every test sample contribute to its own farm mean, which is leakage.

    A farm unseen during fit (possible in a small fold) falls back to the global mean.
    """

    def fit(self, X, y=None):
        Xf, farm = X[:, :-1], X[:, -1]
        self.global_mean_ = Xf.mean(axis=0)
        self.means_ = {f: Xf[farm == f].mean(axis=0) for f in np.unique(farm)}
        return self

    def transform(self, X):
        Xf, farm = X[:, :-1].copy(), X[:, -1]
        for f in np.unique(farm):
            Xf[farm == f] -= self.means_.get(f, self.global_mean_)
        return Xf


def make_pipeline_clf(**overrides):
    """Classifier preceded by within-farm centering; grid keys gain the 'rf__' prefix."""
    return Pipeline([('center', WithinFarmCentering()), ('rf', make_clf(**overrides))])


def make_pipeline_reg(**overrides):
    return Pipeline([('center', WithinFarmCentering()), ('rf', make_reg(**overrides))])


def prefix_grid(grid):
    return {f'rf__{k}': v for k, v in grid.items()}


def attach_farm(X_df):
    """Append the farm code as the last column, for the centering transformer."""
    return np.column_stack([X_df.values, FARM.reindex(X_df.index).values])


def center_y_by_farm(y, farm, means=None):
    """Within-farm centering of a continuous response (fixed-effects 'within' estimator).

    Returns the centered values and the farm means used, so a training fold's means
    can be applied unchanged to its test fold.
    """
    if means is None:
        means = {f: y[farm == f].mean() for f in np.unique(farm)}
    gm = float(np.mean(list(means.values())))
    return np.array([yi - means.get(f, gm) for yi, f in zip(y, farm)]), means


# ── Disk cache: a full run takes hours, a kernel restart must not cost them ──
# CACHE_VERSION is part of the filename. Bump it whenever the procedure changes, so
# results from an earlier definition can never be picked up silently: v1 files came
# from 200 permutations at 500 trees and a differently defined importance step.
CACHE_VERSION = f'v2{OUT_TAG}'


def _cache_path(key, ext='json'):
    return CACHE_DIR / f'{CACHE_VERSION}_{"quick" if QUICK else "full"}_{key}.{ext}'


def cache_load(key, expect=None):
    """Read a cache entry, but only if it was produced with the same settings.

    `expect` names the parameters the caller depends on. A stored entry whose
    '_params' differ is treated as a miss and the reason is printed. Without this
    check, raising n_permutations and re-running silently returns the old p-value --
    which is exactly what happened before the check existed.
    """
    p = _cache_path(key)
    if not p.exists() or FORCE_RECOMPUTE:
        return None
    obj = json.loads(p.read_text())
    if expect:
        stored = obj.get('_params', {})
        diff = {k: (stored.get(k), v) for k, v in expect.items() if stored.get(k) != v}
        if diff:
            msg = ', '.join(f'{k} {was} -> {now}' for k, (was, now) in diff.items())
            print(f'  [stale] {key:44s} {msg}')
            return None
    return obj


def cache_save(key, obj, params=None):
    if params:
        obj = dict(obj, _params=params)
    _cache_path(key).write_text(json.dumps(obj, indent=1))


def imp_to_json(df):
    """Importance DataFrame -> JSON-safe dict (and back via imp_from_json)."""
    return {'features': list(df.index),
            'cols': {c: [float(v) for v in df[c]] for c in df.columns}}


def imp_from_json(obj):
    return pd.DataFrame(obj['cols'], index=obj['features'])


def _estimator_and_grid(kind, blocked, **overrides):
    """Bare model, or the same model behind within-farm centering."""
    if blocked:
        est = make_pipeline_clf(**overrides) if kind == 'clf' else make_pipeline_reg(**overrides)
        return est, prefix_grid(PARAM_GRID)
    est = make_clf(**overrides) if kind == 'clf' else make_reg(**overrides)
    return est, PARAM_GRID


def strip_prefix(params):
    """'rf__min_samples_leaf' -> 'min_samples_leaf', so cached params are comparable."""
    return {k.split('__', 1)[-1]: v for k, v in params.items()}


def nested_cv(X, y, key, kind='clf', blocked=False):
    """Nested CV. Returns mean/std of the outer-fold scores plus the modal
    hyperparameters, i.e. the setting that won most often across outer folds.

    Hyperparameters are selected by the inner GridSearchCV inside every outer fold,
    so the reported score is not inflated by the selection itself. Different outer
    folds may well have chosen different settings; 'best_params' is only the modal
    winner and is used downstream (OOF, permutation test, importances), never for
    the score reported here.
    """
    # repr(), not the structure itself: JSON turns tuples into lists, so a structural
    # comparison would never match after a round-trip and every entry would look stale.
    want = {'n_trees': N_TREES, 'n_repeats': N_REPEATS, 'inner': INNER_SPLITS,
            'grid': repr(sorted(PARAM_GRID.items()))}
    cached = cache_load(key, expect=want)
    if cached is not None:
        print(f'  [cache] {key}')
        return cached

    t0 = time.time()
    est, grid = _estimator_and_grid(kind, blocked)
    if kind == 'clf':
        outer   = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS,
                                          random_state=RANDOM_STATE)
        inner   = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True,
                                  random_state=RANDOM_STATE)
        search  = GridSearchCV(est, grid, cv=inner,
                               scoring='balanced_accuracy', n_jobs=N_JOBS_CV)
        scoring = {'bal_acc': 'balanced_accuracy',
                   'acc':     'accuracy',
                   'mcc':     make_scorer(matthews_corrcoef)}
    else:
        outer   = RepeatedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS,
                                random_state=RANDOM_STATE)
        inner   = KFold(n_splits=INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE)
        search  = GridSearchCV(est, grid, cv=inner, scoring='r2', n_jobs=N_JOBS_CV)
        scoring = {'r2': 'r2', 'neg_mae': 'neg_mean_absolute_error'}

    cv_res = cross_validate(search, X, y, cv=outer, scoring=scoring,
                            return_estimator=True, n_jobs=N_JOBS_CV)

    modal = Counter(tuple(sorted(strip_prefix(e.best_params_).items()))
                    for e in cv_res['estimator'])
    out = {'best_params':   dict(modal.most_common(1)[0][0]),
           'n_outer_folds': int(len(cv_res['test_' + list(scoring)[0]])),
           'runtime_s':     round(time.time() - t0, 1)}
    for m in scoring:
        out[f'{m}_mean'] = float(cv_res[f'test_{m}'].mean())
        out[f'{m}_std']  = float(cv_res[f'test_{m}'].std())
    cache_save(key, out, params=want)
    print(f'  {key:42s} {out["runtime_s"]:6.1f}s  params={out["best_params"]}')
    return out


def oof_predictions(X, y, best_params, key, blocked=False):
    """Single StratifiedKFold pass with fixed hyperparameters: every sample predicted
    exactly once, as a confusion matrix and a calibration curve require. A repeated CV
    predicts each sample several times and cannot serve that purpose.

    Cached so that confusion matrices and calibration curves can be redrawn without
    fitting anything.
    """
    want = {'n_trees': N_TREES, 'best_params': best_params}
    hit = cache_load(key, expect=want)
    if hit is not None:
        return np.array(hit['y_pred']), np.array(hit['y_proba'])

    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    est, _ = _estimator_and_grid('clf', blocked, **best_params)
    y_pred  = cross_val_predict(est, X, y, cv=cv, method='predict')
    y_proba = cross_val_predict(est, X, y, cv=cv, method='predict_proba')
    cache_save(key, {'y_pred': [str(v) for v in y_pred],
                     'y_proba': [[float(p) for p in row] for row in y_proba]}, params=want)
    return y_pred, y_proba


def run_permutation_test(X, y, best_params, key, kind='clf', blocked=False,
                         n_perm=None):
    """Empirical null from shuffled labels -> p-value for the observed score.

    permutation_test_score puts the observed score and every null replicate through
    the SAME procedure, which is what makes the p-value valid. Fewer trees are used
    here than in the main models (N_TREES_PERM) -- allowed for the same reason.
    """
    n_perm = n_perm or N_PERM
    want = {'n_perm': n_perm, 'n_trees': N_TREES_PERM, 'best_params': best_params}
    cached = cache_load(key, expect=want)
    if cached is not None:
        print(f'  [cache] {key}')
        return cached

    t0 = time.time()
    est, _ = _estimator_and_grid(kind, blocked,
                                 n_estimators=N_TREES_PERM, **best_params)
    if kind == 'clf':
        cv, scoring = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                                      random_state=RANDOM_STATE), 'balanced_accuracy'
    else:
        cv, scoring = KFold(n_splits=N_SPLITS, shuffle=True,
                            random_state=RANDOM_STATE), spearman_scorer

    observed, null_scores, pval = permutation_test_score(
        est, X, y, cv=cv, scoring=scoring, n_permutations=n_perm,
        random_state=RANDOM_STATE, n_jobs=N_JOBS_CV)

    out = {'observed': float(observed), 'pvalue': float(pval), 'n_perm': int(n_perm),
           'null': [float(v) for v in null_scores],
           'runtime_s': round(time.time() - t0, 1)}
    cache_save(key, out, params=want)
    print(f'  {key:42s} {out["runtime_s"]:6.1f}s  observed={observed:.3f}  '
          f'p={pval:.4f} (n={n_perm})')
    return out


def fold_importances(X, y, best_params, feat_ids, key, blocked=False):
    """All three importance measures, computed fold-wise instead of on a refit model.

    The previous version fitted one model on all data and measured importance there.
    With min_samples_leaf=1 that model memorises the training set (balanced accuracy
    1.000), so shuffling any single feature changes nothing -- hundreds of correlated
    features offer redundant paths to the same leaf, and permutation importance came
    out exactly zero. Measuring on held-out folds removes that failure mode.

    Returns a DataFrame indexed by feature id -- not by taxonomy label -- with mean
    and std across folds for each measure; the std separates taxa that matter
    consistently from one-fold flukes. Labels are attached when a figure is drawn, so
    re-parsing the taxonomy never costs a recomputation.

    The feature ids are fingerprinted into the cache key: a different feature set (a
    changed prevalence threshold, say) must not be served importances computed for
    other columns just because the hyperparameters happen to match.
    """
    want = {'n_trees': N_TREES, 'n_perm_imp': N_PERM_IMP, 'best_params': best_params,
            'feat': hashlib.md5('|'.join(feat_ids).encode()).hexdigest()[:12]}
    hit = cache_load(key, expect=want)
    if hit is not None:
        return imp_from_json(hit)

    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    gini, shap_v, perm_v = [], [], []

    for tr, te in cv.split(X, y):
        est, _ = _estimator_and_grid('clf', blocked, **best_params)
        est.fit(X[tr], y[tr])
        rf = est.named_steps['rf'] if blocked else est

        gini.append(rf.feature_importances_)

        # SHAP and permutation importance both run on the held-out fold.
        X_te = est.named_steps['center'].transform(X[te]) if blocked else X[te]
        sv = shap.TreeExplainer(rf).shap_values(X_te, check_additivity=False)
        if isinstance(sv, list):                 # shap < 0.42: one array per class
            sv_abs = np.mean([np.abs(s) for s in sv], axis=0)
        elif np.ndim(sv) == 3:                   # shap >= 0.42: (n, features, classes)
            sv_abs = np.abs(sv).mean(axis=2)
        else:
            sv_abs = np.abs(sv)
        shap_v.append(sv_abs.mean(axis=0))

        # Permutation importance runs on the bare forest and the already-centered
        # test matrix, not on the pipeline. Feeding the pipeline would also permute
        # the appended farm column and return one value too many.
        pi = sklearn_perm_imp(rf, X_te, y[te], n_repeats=N_PERM_IMP,
                              random_state=RANDOM_STATE,
                              scoring='balanced_accuracy', n_jobs=N_JOBS_CV)
        perm_v.append(pi.importances_mean)

    out = {}
    for name, vals in [('gini', gini), ('shap', shap_v), ('perm', perm_v)]:
        arr = np.vstack(vals)
        out[f'{name}_mean'] = arr.mean(axis=0)
        out[f'{name}_std']  = arr.std(axis=0)
    df = pd.DataFrame(out, index=feat_ids)
    cache_save(key, imp_to_json(df), params=want)
    return df


def spearman_scorer(estimator, X, y):
    """Rank correlation between prediction and truth. Primary regression metric:
    pi_clsm is zero-inflated and right-skewed, so a rank measure is more honest
    than R^2 here."""
    return spearmanr(y, estimator.predict(X))[0]

# ===== cell 7 ==========================================================
# Which targets get a second, farm-blocked fit.
BLOCKED_TARGETS = ['sample_type', 'pollination', 'fungal_colonization']
VARIANTS = ['raw', 'within-farm']

RESULTS = {}

for target in CLF_TARGETS:
    print(f'\n=== {target} ===')
    variants = VARIANTS if target in BLOCKED_TARGETS else ['raw']
    for variant in variants:
        blocked = variant == 'within-farm'
        for ds in DATASETS:
            X_df, y_s = build_xy(target, ds)
            # Blocked runs carry the farm code as an extra column; the pipeline's
            # first step consumes and strips it.
            X = attach_farm(X_df) if blocked else X_df.values
            y = y_s.values
            key = f'{target}__{ds}' + ('__farmblock' if blocked else '')

            n_perm = (N_PERM_STRONG if (target, ds, variant) in STRONG_COMBOS
                      else N_PERM)

            nested = nested_cv(X, y, key + '__nested', kind='clf', blocked=blocked)
            bp     = nested['best_params']
            y_pred, y_proba = oof_predictions(X, y, bp, key + '__oof', blocked=blocked)
            perm   = run_permutation_test(X, y, bp, key + '__perm', kind='clf',
                                          blocked=blocked, n_perm=n_perm)
            imp    = fold_importances(X, y, bp, FEATURES[ds]['ids'],
                                      key + '__imp', blocked=blocked)

            classes = list(np.unique(y))
            try:
                roc = (roc_auc_score(y, y_proba[:, 1]) if len(classes) == 2
                       else roc_auc_score(y, y_proba, multi_class='ovr',
                                          average='macro', labels=classes))
            except ValueError:
                roc = np.nan

            RESULTS[(target, ds, variant)] = dict(
                nested=nested, perm=perm, imp=imp, classes=classes,
                X_df=X_df, y=y, y_pred=y_pred, y_proba=y_proba,
                oof_bal_acc=balanced_accuracy_score(y, y_pred),
                oof_mcc=matthews_corrcoef(y, y_pred),
                roc_auc=roc, n=len(y), blocked=blocked,
            )
print('\nDone.')

# ===== cell 8 ==========================================================
REG_RESULTS = {}
print('=== pi_clsm (regression) ===')

for variant in VARIANTS:
    blocked = variant == 'within-farm'
    print(f'\n--- {variant} ---')
    for ds in DATASETS:
        X_df, y_s = build_xy('pi_clsm', ds)
        X_raw = X_df.values
        y     = y_s.values.astype(float)
        farm  = FARM.reindex(X_df.index).values
        X     = attach_farm(X_df) if blocked else X_raw
        key   = f'pi_clsm__{ds}' + ('__farmblock' if blocked else '')

        if not blocked:
            nested = nested_cv(X, y, key + '__nested', kind='reg')
            bp     = nested['best_params']
            cv     = KFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
            y_eval = y
            want   = {'n_trees': N_TREES, 'best_params': bp}
            hit    = cache_load(key + '__oofreg', expect=want)
            if hit is not None:
                y_oof = np.array(hit['y_oof'])
            else:
                y_oof = cross_val_predict(make_reg(**bp), X, y, cv=cv)
                cache_save(key + '__oofreg',
                           {'y_oof': [float(v) for v in y_oof],
                            'y_eval': [float(v) for v in y_eval]}, params=want)
            perm   = run_permutation_test(X, y, bp, key + '__perm', kind='reg')
        else:
            # The response is centered as well, giving a proper within (fixed-effects)
            # estimator: "inside one farm, do flowers richer in taxon X carry more
            # pollen?". Farm means of y cannot live in a feature transformer, so the
            # outer loop is manual -- means come from the training fold only, which
            # keeps the centering leak-free exactly as it is on the feature side.
            nested = nested_cv(X, center_y_by_farm(y, farm)[0],
                               key + '__nested', kind='reg', blocked=True)
            bp     = nested['best_params']
            cv     = KFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
            want = {'n_trees': N_TREES, 'best_params': bp}
            hit  = cache_load(key + '__oofreg', expect=want)
            if hit is not None:
                y_oof  = np.array(hit['y_oof'])
                y_eval = np.array(hit['y_eval'])
            else:
                y_eval = np.empty_like(y)
                y_oof  = np.empty_like(y)
                for tr, te in cv.split(X):
                    y_tr_c, means = center_y_by_farm(y[tr], farm[tr])
                    y_te_c, _     = center_y_by_farm(y[te], farm[te], means)
                    model = make_pipeline_reg(**bp).fit(X[tr], y_tr_c)
                    y_oof[te]  = model.predict(X[te])
                    y_eval[te] = y_te_c
                cache_save(key + '__oofreg',
                           {'y_oof': [float(v) for v in y_oof],
                            'y_eval': [float(v) for v in y_eval]}, params=want)
            perm = run_permutation_test(X, center_y_by_farm(y, farm)[0], bp,
                                        key + '__perm', kind='reg', blocked=True)

        rho, _ = spearmanr(y_eval, y_oof)
        REG_RESULTS[(ds, variant)] = dict(
            nested=nested, perm=perm, X_df=X_df, y=y_eval, y_raw=y, y_oof=y_oof,
            spearman=rho, r2=r2_score(y_eval, y_oof),
            mae=mean_absolute_error(y_eval, y_oof), n=len(y), blocked=blocked,
        )
        print(f'  {ds:9s} rho={rho:.3f}  R2={REG_RESULTS[(ds, variant)]["r2"]:.3f}  '
              f'MAE={REG_RESULTS[(ds, variant)]["mae"]:.1f}  p={perm["pvalue"]:.4f}')

# ===== cell 9 ==========================================================
rows = []
for (target, ds, variant), r in RESULTS.items():
    rows.append(dict(
        target=target, dataset=ds, variant=variant, n=r['n'],
        n_classes=len(r['classes']), chance=round(1 / len(r['classes']), 3),
        nested_bal_acc=r['nested']['bal_acc_mean'], nested_std=r['nested']['bal_acc_std'],
        nested_mcc=r['nested']['mcc_mean'],
        oof_bal_acc=r['oof_bal_acc'], oof_mcc=r['oof_mcc'], roc_auc=r['roc_auc'],
        perm_observed=r['perm']['observed'], perm_p=r['perm']['pvalue'],
        n_perm=r['perm'].get('n_perm', np.nan),
        best_params=str(r['nested']['best_params']),
    ))
summary_df = pd.DataFrame(rows)

# Balanced accuracy expressed as progress from chance to perfect. Targets have
# different chance levels (0.143 / 0.333 / 0.5), so raw accuracy is not comparable
# across rows and a single colour scale over it would mislead.
summary_df['above_chance'] = ((summary_df['nested_bal_acc'] - summary_df['chance'])
                              / (1 - summary_df['chance']))

for (ds, variant), r in REG_RESULTS.items():
    summary_df.loc[len(summary_df)] = dict(
        target='pi_clsm (reg)', dataset=ds, variant=variant, n=r['n'],
        n_classes=np.nan, chance=np.nan, nested_bal_acc=np.nan, nested_std=np.nan,
        nested_mcc=np.nan, oof_bal_acc=np.nan, oof_mcc=np.nan, roc_auc=np.nan,
        perm_observed=r['perm']['observed'], perm_p=r['perm']['pvalue'],
        n_perm=r['perm'].get('n_perm', np.nan), above_chance=np.nan,
        best_params=str(r['nested']['best_params']))

summary_df.to_csv(TAB_DIR / f'{FIG_PREFIX}rf_summary.csv', index=False)

show = summary_df.copy()
for c in ['nested_bal_acc', 'nested_std', 'nested_mcc', 'roc_auc',
          'perm_observed', 'above_chance']:
    show[c] = show[c].map(lambda v: '—' if pd.isna(v) else f'{v:.3f}')
show['perm_p'] = show['perm_p'].map(lambda v: f'{v:.4f}')
print(show[['target', 'dataset', 'variant', 'n', 'chance', 'nested_bal_acc',
            'nested_std', 'above_chance', 'nested_mcc', 'roc_auc',
            'perm_p', 'n_perm']].to_string(index=False))

print('\nRegression (pi_clsm):')
for (ds, variant), r in REG_RESULTS.items():
    print(f'  {ds:9s} {variant:12s} n={r["n"]}  Spearman={r["spearman"]:.3f}  '
          f'R2={r["r2"]:.3f}  MAE={r["mae"]:.1f}  p={r["perm"]["pvalue"]:.4f}')

# ===== cell 10 ==========================================================
# ── Heatmap: how far each model gets from chance towards perfect ────────────
# Colour encodes (bal_acc - chance) / (1 - chance): 0 = chance, 1 = perfect. The cell
# annotation keeps the raw balanced accuracy so nothing is hidden behind the scaling.
clf_rows = summary_df[summary_df['target'] != 'pi_clsm (reg)'].copy()
clf_rows['row'] = clf_rows['target'] + np.where(
    clf_rows['variant'] == 'within-farm', '  (within-farm)', '')
row_order = [f'{t}{s}' for t in CLF_TARGETS
             for s in ['', '  (within-farm)']
             if f'{t}{s}' in set(clf_rows['row'])]

pivot = clf_rows.pivot(index='row', columns='dataset',
                       values='above_chance').reindex(index=row_order, columns=DATASETS)
annot = clf_rows.pivot(index='row', columns='dataset',
                       values='nested_bal_acc').reindex(index=row_order, columns=DATASETS)

fig, ax = plt.subplots(figsize=(7, 0.55 * len(row_order) + 2))
sns.heatmap(pivot.astype(float), annot=annot.map(lambda v: f'{v:.3f}'), fmt='',
            cmap='RdYlGn', vmin=0, vmax=1,          # 0 = chance, 1 = perfect
            linewidths=0.5, ax=ax,
            cbar_kws={'label': '(bal. acc - chance) / (1 - chance)'})
ax.set_title('Nested-CV performance above chance\n'
             'cell annotation = raw balanced accuracy', fontweight='bold', fontsize=11)
ax.set_xlabel(''); ax.set_ylabel('')
plt.tight_layout()
fig.savefig(FIG_DIR / f'{FIG_PREFIX}summary_heatmap.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()

# ===== cell 11 ==========================================================
def plot_target_panel(target, variant='raw', n_top=N_TOP_IMP):
    """Two rows: confusion matrices on top, top-N Gini importances below.

    height_ratios gives the bar row more space than the (square, hence small)
    confusion matrices; tick label size drops to 7 because taxon labels are long.
    Error bars on the bars are the spread across CV folds -- they separate taxa that
    matter consistently from ones that surfaced in a single fold.
    """
    order = TARGETS[target]['order']
    fig, axes = plt.subplots(2, len(DATASETS), figsize=(5.4 * len(DATASETS), 8.8),
                             gridspec_kw={'height_ratios': [1, 1.7]})
    suffix = '' if variant == 'raw' else '  —  within-farm (between-farm variation removed)'
    fig.suptitle(f'{target}{suffix}', fontsize=13, fontweight='bold', y=0.995)

    for col, ds in enumerate(DATASETS):
        r  = RESULTS[(target, ds, variant)]
        cn = order if order is not None else r['classes']
        cm = confusion_matrix(r['y'], r['y_pred'], labels=cn)
        ConfusionMatrixDisplay(cm, display_labels=cn).plot(
            ax=axes[0, col], colorbar=False, cmap='Blues')
        roc = '—' if np.isnan(r['roc_auc']) else f'{r["roc_auc"]:.3f}'
        axes[0, col].set_title(
            f'{ds}  (n={r["n"]})\n'
            f'nested bal.acc = {r["nested"]["bal_acc_mean"]:.3f} '
            f'± {r["nested"]["bal_acc_std"]:.3f}\n'
            f'MCC = {r["nested"]["mcc_mean"]:.3f}   ROC-AUC = {roc}\n'
            f'permutation p = {r["perm"]["pvalue"]:.4f} '
            f'(n={r["perm"].get("n_perm", N_PERM)})', fontsize=9)
        axes[0, col].tick_params(axis='x', rotation=30, labelsize=8)

        # Select and align on feature ids, then rename to taxonomy labels for the
        # axis -- doing it the other way round would break the reindex.
        imp = r['imp']['gini_mean'].nlargest(n_top).sort_values()
        err = r['imp']['gini_std'].reindex(imp.index)
        imp = imp.rename(index=FEATURES[ds]['label'])
        axes[1, col].barh(imp.index, imp.values, xerr=err.values,
                          color=DS_COLOURS[ds], alpha=0.85,
                          error_kw=dict(lw=0.6, ecolor='grey'))
        axes[1, col].set_title(f'{ds} — top {n_top} Gini importance', fontsize=9)
        axes[1, col].set_xlabel('Mean decrease in impurity (mean ± SD over folds)')
        axes[1, col].tick_params(axis='y', labelsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    tag = '' if variant == 'raw' else '_farmblock'
    fig.savefig(FIG_DIR / f'{FIG_PREFIX}panel_{target}{tag}.png', dpi=FIG_DPI,
                bbox_inches='tight', facecolor='white')
    plt.show()

for t in CLF_TARGETS:
    for v in (VARIANTS if t in BLOCKED_TARGETS else ['raw']):
        plot_target_panel(t, v)

# ===== cell 12 ==========================================================
perm_rows = [(t, v) for t in CLF_TARGETS
             for v in (VARIANTS if t in BLOCKED_TARGETS else ['raw'])]

fig, axes = plt.subplots(len(perm_rows), len(DATASETS),
                         figsize=(4.6 * len(DATASETS), 3.1 * len(perm_rows)))
fig.suptitle('Permutation null distributions', fontsize=12, fontweight='bold', y=0.995)

for row, (target, variant) in enumerate(perm_rows):
    for col, ds in enumerate(DATASETS):
        ax   = axes[row, col]
        r    = RESULTS[(target, ds, variant)]
        null = np.asarray(r['perm']['null'])

        ax.hist(null, bins=25, color='grey', alpha=0.6, edgecolor='white')
        # Dashed vertical = observed score; dotted = theoretical chance level.
        ax.axvline(r['perm']['observed'], color=DS_COLOURS[ds], lw=2, ls='--')
        ax.axvline(1 / len(r['classes']), color='black', lw=1, ls=':')
        tag = '' if variant == 'raw' else ' | within-farm'
        ax.set_title(f'{target}{tag} | {ds}\n'
                     f'observed = {r["perm"]["observed"]:.3f}, '
                     f'p = {r["perm"]["pvalue"]:.4f} '
                     f'({r["perm"].get("n_perm", N_PERM)} shuffles)', fontsize=8.5)
        ax.set_xlabel('Balanced accuracy'); ax.set_ylabel('Permutations')
        ax.tick_params(labelsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(FIG_DIR / f'{FIG_PREFIX}permutation_tests.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()

# ===== cell 13 ==========================================================
BINARY_TARGETS = [t for t in CLF_TARGETS if len(TARGETS[t]['order'] or []) == 2]

fig, axes = plt.subplots(1, len(BINARY_TARGETS),
                         figsize=(5.2 * len(BINARY_TARGETS), 4.6), squeeze=False)
fig.suptitle('Calibration of out-of-fold probabilities', fontsize=12, fontweight='bold')

for col, target in enumerate(BINARY_TARGETS):
    ax = axes[0, col]
    ax.plot([0, 1], [0, 1], 'k:', lw=1, label='perfect')   # reference diagonal
    for ds in DATASETS:
        r = RESULTS[(target, ds, 'raw')]
        pos = r['classes'][1]                       # probability of the second class
        y_bin = (r['y'] == pos).astype(int)
        p_hat = r['y_proba'][:, 1]
        # n_bins=8 keeps ~8-13 samples per bin at these sample sizes.
        frac, mean_pred = calibration_curve(y_bin, p_hat, n_bins=8, strategy='quantile')
        ax.plot(mean_pred, frac, 'o-', color=DS_COLOURS[ds], alpha=0.85,
                label=f'{ds} (Brier={brier_score_loss(y_bin, p_hat):.3f})')
    ax.set_title(f'{target}\npositive class = '
                 f'"{RESULTS[(target, DATASETS[0], "raw")]["classes"][1]}"', fontsize=10)
    ax.set_xlabel('Mean predicted probability'); ax.set_ylabel('Observed frequency')
    ax.legend(fontsize=8, loc='upper left')

plt.tight_layout()
fig.savefig(FIG_DIR / f'{FIG_PREFIX}calibration.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()

# ===== cell 14 ==========================================================
POLL_COLOURS = {'none': '#4C72B0', 'weak': '#DD8452', 'strong': '#55A868'}

# One row per variant. pi_clsm spans 0 to ~350 with most flowers near zero, so the
# raw row is shown on a log1p scale where the bulk is actually readable. The
# within-farm row is already centered (values straddle zero) and stays linear.
fig, axes = plt.subplots(len(VARIANTS), len(DATASETS),
                         figsize=(5.2 * len(DATASETS), 4.6 * len(VARIANTS)))
fig.suptitle('pi_clsm — out-of-fold RF regression', fontsize=12,
             fontweight='bold', y=0.995)

for row, variant in enumerate(VARIANTS):
    for col, ds in enumerate(DATASETS):
        ax  = axes[row, col]
        r   = REG_RESULTS[(ds, variant)]
        # Colour by the 3-class grouping of the RAW pollen load, so the two rows
        # stay visually comparable even though the y axis differs.
        grp = pd.Series(r['y_raw']).apply(
            lambda v: 'none' if v == 0 else ('weak' if v < PI_THRESH_POLL else 'strong'))

        if variant == 'raw':
            xs, ys = np.log1p(r['y']), np.log1p(np.clip(r['y_oof'], 0, None))
            xlab = 'Observed pi_clsm (log1p grains/stigma)'
        else:
            xs, ys = r['y'], r['y_oof']
            xlab = 'Observed pi_clsm, farm-centered'

        for g, c in POLL_COLOURS.items():
            m = (grp == g).values
            ax.scatter(xs[m], ys[m], c=c, s=36, alpha=0.8,
                       edgecolors='black', linewidths=0.3, label=g)
        lim = [min(xs.min(), ys.min()), max(xs.max(), ys.max())]
        ax.plot(lim, lim, 'k:', lw=1)                  # 1:1 line, not a fit
        if variant == 'within-farm':
            ax.axhline(0, color='grey', lw=0.6); ax.axvline(0, color='grey', lw=0.6)

        ax.set_title(f'{ds} | {variant}  (n={r["n"]})\n'
                     f'Spearman rho = {r["spearman"]:.3f}  '
                     f'(p = {r["perm"]["pvalue"]:.4f})\n'
                     f'R² = {r["r2"]:.3f}   MAE = {r["mae"]:.1f}', fontsize=9)
        ax.set_xlabel(xlab)
        if col == 0:
            ax.set_ylabel('Predicted (out-of-fold)')
        if row == 0 and col == 0:
            ax.legend(title='3-class group', fontsize=8, title_fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(FIG_DIR / f'{FIG_PREFIX}regression_pi_clsm.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()

# ===== cell 15 ==========================================================
# ── Diagnostic: how strongly is library size tied to each target? ───────────
lib_rows = []
for ds, counts_raw in [('ITS1', counts_its1), ('16S', counts_16s)]:
    depth = counts_raw.sum(axis=0)                      # reads per sample
    for target in CLF_TARGETS:
        _, y_s = build_xy(target, ds)
        d = depth.reindex(y_s.index)
        # eta^2 from a one-way layout: fraction of depth variance explained by the class
        grand = d.mean()
        ss_between = sum(len(g) * (g.mean() - grand) ** 2 for _, g in d.groupby(y_s.values))
        ss_total   = ((d - grand) ** 2).sum()
        lib_rows.append(dict(dataset=ds, target=target,
                             eta2_depth=ss_between / ss_total,
                             median_depth=int(d.median())))
lib_df = pd.DataFrame(lib_rows)
print('Variance in sequencing depth explained by each target (eta^2):')
print(lib_df.to_string(index=False))
print('\nHigh eta^2 means raw counts give the forest a technical shortcut.')

# ===== cell 16 ==========================================================
TRANSFORMS = ['raw', 'tss', 'clr', 'rclr']

def score_fixed(X_df, y_s):
    """Single 5-fold CV with default hyperparameters — used by the sensitivity panels."""
    if y_s.nunique() < 2 or len(y_s) < 2 * N_SPLITS:
        return np.nan
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    pred = cross_val_predict(make_clf(), X_df.values, y_s.values, cv=cv)
    return balanced_accuracy_score(y_s.values, pred)

tr_rows = []
for ds, counts_raw, tax_d in [('ITS1', counts_its1, tax_its1), ('16S', counts_16s, tax_16s)]:
    for how in TRANSFORMS:
        X_t, _ = build_features(counts_raw, meta, tax_d, how=how)
        for target in CLF_TARGETS:
            X_sub, y_sub = build_xy(target, ds, X_override=X_t)
            tr_rows.append(dict(dataset=ds, transform=how, target=target,
                                bal_acc=score_fixed(X_sub, y_sub)))
        print(f'  {ds} {how} done')

tr_df = pd.DataFrame(tr_rows)
tr_df.to_csv(TAB_DIR / f'{FIG_PREFIX}transform_comparison.csv', index=False)

# ===== cell 17 ==========================================================
fig, axes = plt.subplots(2, len(CLF_TARGETS), figsize=(4.1 * len(CLF_TARGETS), 7.4),
                         sharey=True)
fig.suptitle('Sensitivity to the feature transformation (fixed hyperparameters, 5-fold CV)',
             fontsize=12, fontweight='bold')

# Grey bars for the alternatives, the dataset colour for CLR (the default in use).
for row, ds in enumerate(['ITS1', '16S']):
    for col, target in enumerate(CLF_TARGETS):
        ax  = axes[row, col]
        sub = tr_df[(tr_df.dataset == ds) & (tr_df.target == target)].set_index('transform')
        sub = sub.reindex(TRANSFORMS)
        colours = [DS_COLOURS[ds] if t == 'clr' else '#BBBBBB' for t in TRANSFORMS]
        ax.bar(range(len(TRANSFORMS)), sub['bal_acc'], color=colours, alpha=0.9)
        n_cls = len(TARGETS[target]['order'] or RESULTS[(target, ds, 'raw')]['classes'])
        ax.axhline(1 / n_cls, color='black', ls=':', lw=1)     # chance level
        ax.set_xticks(range(len(TRANSFORMS)))
        ax.set_xticklabels(TRANSFORMS, fontsize=9)
        ax.set_title(f'{ds} | {target}', fontsize=9, fontweight='bold')
        if col == 0:
            ax.set_ylabel('Balanced accuracy')
        ax.set_ylim(0, 1.02)

plt.tight_layout()
fig.savefig(FIG_DIR / f'{FIG_PREFIX}transform_comparison.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()

# ===== cell 18 ==========================================================
unbagged_ids = list(meta.index[meta['sample_type'] == 'unbagged_flower'])
PREV_CONFIGS = [('0 %', 0.00, None), ('5 %', 0.05, None), ('10 % (used)', 0.10, None),
                ('20 %', 0.20, None), ('10 % within unbagged', 0.10, unbagged_ids)]

pv_rows = []
for ds, counts_raw, tax_d in [('ITS1', counts_its1, tax_its1), ('16S', counts_16s, tax_16s)]:
    for label, thr, subset in PREV_CONFIGS:
        X_t, _ = build_features(counts_raw, meta, tax_d, how='clr',
                                thresh=thr, subset_ids=subset)
        for target in CLF_TARGETS:
            X_sub, y_sub = build_xy(target, ds, X_override=X_t)
            pv_rows.append(dict(dataset=ds, filter=label, target=target,
                                n_features=X_t.shape[1], bal_acc=score_fixed(X_sub, y_sub)))
        print(f'  {ds} {label}: {X_t.shape[1]} features')

pv_df = pd.DataFrame(pv_rows)
pv_df.to_csv(TAB_DIR / f'{FIG_PREFIX}prevalence_comparison.csv', index=False)

# ===== cell 19 ==========================================================
fig, axes = plt.subplots(2, len(CLF_TARGETS), figsize=(4.3 * len(CLF_TARGETS), 7.8),
                         sharey=True)
fig.suptitle('Effect of the prevalence threshold (fixed hyperparameters, 5-fold CV)',
             fontsize=12, fontweight='bold')
labels = [c[0] for c in PREV_CONFIGS]

for row, ds in enumerate(['ITS1', '16S']):
    for col, target in enumerate(CLF_TARGETS):
        ax  = axes[row, col]
        sub = (pv_df[(pv_df.dataset == ds) & (pv_df.target == target)]
               .set_index('filter').reindex(labels))
        # The threshold actually in use is highlighted; the rest are grey.
        colours = [DS_COLOURS[ds] if l.startswith('10 % (') else '#BBBBBB' for l in labels]
        ax.bar(range(len(labels)), sub['bal_acc'], color=colours, alpha=0.9)
        n_cls = len(TARGETS[target]['order'] or RESULTS[(target, ds, 'raw')]['classes'])
        ax.axhline(1 / n_cls, color='black', ls=':', lw=1)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=7.5)
        ax.set_title(f'{ds} | {target}', fontsize=9, fontweight='bold')
        if col == 0:
            ax.set_ylabel('Balanced accuracy')
        ax.set_ylim(0, 1.02)
        # Secondary axis: number of features surviving each threshold.
        ax2 = ax.twinx()
        ax2.plot(range(len(labels)), sub['n_features'], 'ko--', ms=4, lw=1)
        ax2.set_ylabel('# features', fontsize=7, color='grey')
        ax2.tick_params(axis='y', labelsize=7, colors='grey')

plt.tight_layout()
fig.savefig(FIG_DIR / f'{FIG_PREFIX}prevalence_comparison.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()

# ===== cell 20 ==========================================================
def plot_importance_panel(measure, title, xlabel, fname):
    """One row per (target, variant), one column per feature set."""
    rows = [(t, v) for t in CLF_TARGETS
            for v in (VARIANTS if t in BLOCKED_TARGETS else ['raw'])]
    fig, axes = plt.subplots(len(rows), len(DATASETS),
                             figsize=(5.4 * len(DATASETS), 5.0 * len(rows)))
    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.995)

    for row, (target, variant) in enumerate(rows):
        for col, ds in enumerate(DATASETS):
            ax  = axes[row, col]
            imp = RESULTS[(target, ds, variant)]['imp']
            top = imp[f'{measure}_mean'].nlargest(N_TOP_IMP).sort_values()
            err = imp[f'{measure}_std'].reindex(top.index)   # aligned on ids
            top = top.rename(index=FEATURES[ds]['label'])    # labelled for the axis
            ax.barh(top.index, top.values, xerr=err.values,
                    color=DS_COLOURS[ds], alpha=0.85,
                    error_kw=dict(lw=0.6, ecolor='grey'))
            ax.axvline(0, color='grey', lw=0.6)
            tag = '' if variant == 'raw' else ' | within-farm'
            ax.set_title(f'{target}{tag} | {ds}', fontsize=9, fontweight='bold')
            ax.set_xlabel(xlabel)
            ax.tick_params(axis='y', labelsize=7)

    # rect leaves room for suptitle, which tight_layout alone would let the top row cover
    plt.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(FIG_DIR / fname, dpi=FIG_DPI, bbox_inches='tight', facecolor='white')
    plt.show()


plot_importance_panel(
    'perm', f'Permutation importance, out-of-fold (n_repeats={N_PERM_IMP})',
    'Drop in balanced accuracy (mean ± SD over folds)',
    f'{FIG_PREFIX}permutation_importance.png')

plot_importance_panel(
    'shap', 'SHAP importance, out-of-fold (mean |SHAP| over classes)',
    'Mean |SHAP value| (mean ± SD over folds)',
    f'{FIG_PREFIX}shap_importance.png')

# Long-form table of every measure, for cross-referencing with ANCOMBC2 / ISA hits.
# 'feature_id' is the join key -- the same identifier those scripts report -- and
# 'feature_label' is the human-readable name used on the figures.
imp_long = pd.concat(
    {k: v for k, v in
     ((f'{t}|{d}|{v}', RESULTS[(t, d, v)]['imp']) for (t, d, v) in RESULTS)},
    names=['target_dataset_variant', 'feature_id']).reset_index()
imp_long.insert(2, 'feature_label', [
    FEATURES[combo.split('|')[1]]['label'].get(fid, fid)
    for combo, fid in zip(imp_long['target_dataset_variant'], imp_long['feature_id'])])
imp_long.to_csv(TAB_DIR / f'{FIG_PREFIX}importances.csv', index=False)
print('Saved importances table.')

# ===== cell 21 ==========================================================
def imp_by_rank(imp_series, rank_map):
    """Sum fold-averaged Gini importance per deepest-assigned taxonomic rank.

    imp_series is indexed by feature id and rank_map is keyed by feature id, so the
    two are joined on the index. An earlier version reattached the ids positionally,
    which would have silently mismatched had the feature set ever changed.
    """
    key = imp_series.index.map(lambda i: rank_map.get(i, str(i)[:12]))
    return imp_series.groupby(key).sum().sort_values(ascending=False)


def rank_heatmap(ds, n_top=25):
    cols = {}
    for target in CLF_TARGETS:
        r = RESULTS[(target, ds, 'raw')]
        cols[target] = imp_by_rank(r['imp']['gini_mean'], FEATURES[ds]['rank'])
    df  = pd.DataFrame(cols).fillna(0)
    df  = df.loc[df.max(axis=1).nlargest(n_top).index]
    dfn = df.div(df.max(axis=0).replace(0, 1))      # per-column normalisation

    fig, ax = plt.subplots(figsize=(1.9 * len(CLF_TARGETS) + 3, max(6, n_top * 0.36)))
    sns.heatmap(dfn, annot=df.map(lambda v: f'{v:.3f}'), fmt='',
                cmap='YlOrRd', linewidths=0.3, ax=ax,
                cbar_kws={'label': 'Gini importance, normalised per target'})
    ax.set_title(f'{ds} — importance by deepest assigned rank', fontweight='bold')
    ax.set_xlabel(''); ax.set_ylabel('')
    plt.setp(ax.get_xticklabels(), rotation=25, ha='right', fontsize=8)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f'{FIG_PREFIX}rank_importance_{ds}.png', dpi=FIG_DPI,
                bbox_inches='tight', facecolor='white')
    plt.show()
    return df

# Combined has no rank map of its own (its columns carry marker prefixes), so only
# the single-marker feature sets are shown.
for ds in ['ITS1', '16S']:
    rank_heatmap(ds).to_csv(TAB_DIR / f'{FIG_PREFIX}rank_importance_{ds}.csv')
