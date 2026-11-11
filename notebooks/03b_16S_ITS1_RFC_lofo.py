#!/usr/bin/env python
# coding: utf-8
"""03b - Random Forest, leave-one-farm-out (LOFO): bagged vs unbagged.

Can a model trained on six farms tell bagged from unbagged flowers on the seventh,
whose trees and extraction batches it has never seen?

Why this exists. The RFC of notebook 03 uses stratified folds over flowers, so
flowers of the test farm (and of the same tree) are always in the training set
too. Within a farm, bagged and unbagged flowers were extracted in separate
batches; a model can therefore learn the batch of a known farm. Here the test
farm is held out completely. A random batch effect cannot carry over to a farm
the model never saw; a systematic order effect (bagged always extracted first)
still could.

Kept apart from notebook 03 on purpose: own script, own output folder
(results/<version>/rfc_lofo/), own file prefix (03b_lofo_). Nothing here reads or
writes the notebook 03 outputs, except the old nested-CV score drawn as a
reference line in the figure.

Design
  Target    sample_type (bagged / unbagged), all flowers of the sample set
  Features  ITS1 OTU97, 16S ASV, Combined; CLR matrices of 08c (10% prevalence
            on all flowers, CZM zero replacement, as in 08a/08b), exactly as in
            notebook 03
  Variants  raw           the features as they are
            farm-centred  each farm's mean subtracted from its own flowers before
                          the model sees them; uses features only, never labels,
                          so it removes the farm offset without leaking the target
  Outer     LeaveOneGroupOut over farms: 7 folds, one farm held out each time
  Inner     GridSearchCV with LeaveOneGroupOut over the 6 training farms, so
            hyperparameters are also chosen on unseen farms
  Score     balanced accuracy and ROC-AUC on the held-out farm; reported per farm
            and pooled over all held-out predictions
  p-value   permutation_test_score with the same LOFO splits; treatment labels
            are shuffled within farm (groups = farm), modal hyperparameters
  Saved     every flower's held-out prediction, for later checks (for example the
            extraction batch 2024-09-20)

Environment
  CFM_VERSION     v1 | v2 (default v2)
  CFM_SAMPLE_SET  main (288, default) | with_deviations (294)
  LOFO_QUICK      1 = quick pass to test the code (small grid, few trees, 20 shuffles)

Run from the notebooks/ folder (relative paths), on Hamilton via
scripts/03b_hpc_rfc_lofo.sh.
"""

import os
import time
import warnings
from collections import Counter
from pathlib import Path

import biom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut, permutation_test_score

from utils_rfc import load_clr   # CLR matrices of scripts/08c_16S_ITS1_clr_export.R

warnings.filterwarnings('ignore')

# ── 0. Parameters ────────────────────────────────────────────────────────────
CFM_VERSION = os.environ.get('CFM_VERSION', 'v2')
SAMPLE_SET  = os.environ.get('CFM_SAMPLE_SET', 'main')
QUICK       = os.environ.get('LOFO_QUICK', '0') == '1'
assert CFM_VERSION in ('v1', 'v2') and SAMPLE_SET in ('main', 'with_deviations')

ROOT        = Path('..').resolve()
VER_DIR     = ROOT / 'results' / CFM_VERSION
BIOM_DIR    = VER_DIR / 'biom'
BIOM_SUFFIX = '' if SAMPLE_SET == 'main' else '_with_deviations'
SET_TAG     = '_288' if SAMPLE_SET == 'main' else '_294'
OUT_DIR     = VER_DIR / 'rfc_lofo'                 # never the notebook 03 folders
TAB_DIR     = OUT_DIR / 'tables'
FIG_DIR     = OUT_DIR / 'figures'
PREFIX      = '03b_lofo_' + ('quick_' if QUICK else '')
for d in (TAB_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
CLR_CONFIG   = 'all10'       # 08c: 10% prevalence on all flowers, as notebook 03
N_JOBS       = -1            # parallel loops, single-threaded forests (as notebook 03)
FIG_DPI      = 300

if QUICK:
    N_TREES, N_TREES_PERM, N_PERM = 100, 50, 20
    PARAM_GRID = {'min_samples_leaf': [1, 5], 'max_depth': [None]}
else:
    N_TREES, N_TREES_PERM, N_PERM = 500, 200, 1000   # as notebook 03 (N_PERM_STRONG)
    PARAM_GRID = {'min_samples_leaf': [1, 2, 5, 10], 'max_depth': [None, 10, 20]}

CLASSES   = ['bagged_flower', 'unbagged_flower']
POSITIVE  = 'unbagged_flower'
DATASETS  = ['ITS1', '16S', 'Combined']
VARIANTS  = ['raw', 'farm-centred']

print(f'Profile {"QUICK" if QUICK else "FULL"} | data {CFM_VERSION}, set {SAMPLE_SET} '
      f'| trees {N_TREES}, perm {N_PERM} @ {N_TREES_PERM} trees -> {OUT_DIR}')

# ── 1. Data (same feature construction as notebook 03) ───────────────────────
def load_counts(stem):
    """features x samples count table from the biom export of 01a/01b."""
    return biom.load_table(str(BIOM_DIR / f'{stem}{BIOM_SUFFIX}_feature-table.biom')) \
               .to_dataframe(dense=True)


# CLR features from R (08c); the biom counts only serve the consistency check
X_its1 = load_clr('ITS1', VER_DIR, SAMPLE_SET, CLR_CONFIG,
                  counts=load_counts('cfm_its1_otu97_fungi_biosamples'))
X_16s  = load_clr('16S', VER_DIR, SAMPLE_SET, CLR_CONFIG,
                  counts=load_counts('cfm_16s_bacteria_biosamples'))
common = X_its1.index.intersection(X_16s.index)
X_all  = {
    'ITS1':     X_its1,
    '16S':      X_16s,
    'Combined': pd.concat([X_its1.loc[common].add_prefix('ITS1_'),
                           X_16s.loc[common].add_prefix('S16_')], axis=1),
}

meta = pd.read_csv(ROOT / 'data' / 'sample_metadata_with_microscopy.txt',
                   sep='\t', index_col=0, dtype=str)
meta = meta.drop(index='#q2:types', errors='ignore')
meta = meta[meta['sample_type'].isin(CLASSES)]
meta = meta.loc[meta.index.intersection(X_its1.index.union(X_16s.index))]

for ds, X in X_all.items():
    print(f'{ds:9s} {X.shape[1]:4d} features x {X.shape[0]:3d} flowers')


def farm_centre(X, farm):
    """Subtract each farm's own feature mean; labels are not used."""
    return X - X.groupby(farm.values).transform('mean')


# ── 2. LOFO with farm-wise inner tuning ──────────────────────────────────────
def make_clf(**params):
    return RandomForestClassifier(n_estimators=N_TREES, max_features='sqrt',
                                  class_weight='balanced', random_state=RANDOM_STATE,
                                  n_jobs=1, **params)


def run_lofo(X, y, farm):
    """Outer: one farm held out. Inner: grid search, also leave-one-farm-out."""
    rows, preds = [], []
    for train, test in LeaveOneGroupOut().split(X, y, groups=farm):
        search = GridSearchCV(make_clf(), PARAM_GRID, cv=LeaveOneGroupOut(),
                              scoring='balanced_accuracy', n_jobs=N_JOBS)
        search.fit(X[train], y[train], groups=farm[train])
        proba  = search.predict_proba(X[test])[:, list(search.classes_).index(POSITIVE)]
        pred   = search.predict(X[test])
        held   = farm[test][0]
        rows.append(dict(farm=held, n=len(test),
                         bal_acc=balanced_accuracy_score(y[test], pred),
                         roc_auc=roc_auc_score(y[test] == POSITIVE, proba),
                         params=dict(search.best_params_)))
        preds.append(pd.DataFrame({'idx': test, 'pred': pred, 'proba_unbagged': proba}))
    return pd.DataFrame(rows), pd.concat(preds).sort_values('idx')


results, per_farm, predictions = [], [], []

for variant in VARIANTS:
    for ds in DATASETS:
        t0     = time.time()
        X_df   = X_all[ds]
        m      = meta.loc[X_df.index]
        farm   = m['farm_id'].values
        y      = m['sample_type'].values
        if variant == 'farm-centred':
            X_df = farm_centre(X_df, m['farm_id'])
        X = X_df.values

        farms_df, pred_df = run_lofo(X, y, farm)

        # Modal hyperparameters across the 7 outer folds, for the permutation test
        modal = Counter(tuple(sorted(p.items())) for p in farms_df['params'])
        best  = dict(modal.most_common(1)[0][0])

        clf_perm = RandomForestClassifier(n_estimators=N_TREES_PERM, max_features='sqrt',
                                          class_weight='balanced',
                                          random_state=RANDOM_STATE, n_jobs=1, **best)
        observed, null, pval = permutation_test_score(
            clf_perm, X, y, groups=farm, cv=LeaveOneGroupOut(),
            scoring='balanced_accuracy', n_permutations=N_PERM,
            random_state=RANDOM_STATE, n_jobs=N_JOBS)

        y_pred = pred_df['pred'].values
        proba  = pred_df['proba_unbagged'].values
        results.append(dict(
            dataset=ds, variant=variant, n=len(y), n_features=X.shape[1],
            pooled_bal_acc=balanced_accuracy_score(y, y_pred),
            pooled_mcc=matthews_corrcoef(y, y_pred),
            pooled_roc_auc=roc_auc_score(y == POSITIVE, proba),
            mean_farm_bal_acc=farms_df['bal_acc'].mean(),
            min_farm_bal_acc=farms_df['bal_acc'].min(),
            farms_above_chance=int((farms_df['bal_acc'] > 0.5).sum()),
            perm_observed=observed, perm_p=pval, n_perm=N_PERM,
            perm_null_mean=float(np.mean(null)), best_params=str(best),
            runtime_min=round((time.time() - t0) / 60, 1)))
        per_farm.append(farms_df.assign(dataset=ds, variant=variant,
                                        params=farms_df['params'].astype(str)))
        predictions.append(pd.DataFrame({
            'sample_id': X_df.index, 'farm_id': farm, 'tree_id': m['tree_id'].values,
            'extraction_batch': m['extraction_batch'].values, 'true': y,
            'pred': y_pred, 'proba_unbagged': proba, 'dataset': ds, 'variant': variant}))
        r = results[-1]
        print(f'{variant:12s} {ds:9s} pooled bal_acc {r["pooled_bal_acc"]:.3f}  '
              f'AUC {r["pooled_roc_auc"]:.3f}  farms > 0.5: {r["farms_above_chance"]}/7  '
              f'perm p {pval:.4f}  ({r["runtime_min"]} min)', flush=True)

summary = pd.DataFrame(results)
per_farm_df = pd.concat(per_farm)[['dataset', 'variant', 'farm', 'n', 'bal_acc', 'roc_auc', 'params']]
pred_all = pd.concat(predictions)

summary.to_csv(TAB_DIR / f'{PREFIX}summary{SET_TAG}.csv', index=False)
per_farm_df.to_csv(TAB_DIR / f'{PREFIX}per_farm{SET_TAG}.csv', index=False)
pred_all.to_csv(TAB_DIR / f'{PREFIX}predictions{SET_TAG}.csv', index=False)

# ── 3. Figure: balanced accuracy on each held-out farm ───────────────────────
# Reference line: the old notebook 03 nested-CV score (random folds over flowers),
# read only, to show how much of it survives an unseen farm.
old_path = VER_DIR / ('tables' if SAMPLE_SET == 'main' else 'rfc_with_deviations/tables') / '03_hpc_rf_summary.csv'
old = {}
if old_path.exists():
    o = pd.read_csv(old_path)
    o = o[(o['target'] == 'sample_type') & (o['variant'] == 'raw')]
    old = dict(zip(o['dataset'], o['nested_bal_acc']))

DS_COLOURS = {'ITS1': '#9467BD', '16S': '#D62728', 'Combined': '#2CA02C'}  # as notebook 03
fig, axes = plt.subplots(1, len(VARIANTS), figsize=(11, 4.5), sharey=True)   # panel size
for ax, variant in zip(axes, VARIANTS):
    for i, ds in enumerate(DATASETS):
        d = per_farm_df[(per_farm_df['variant'] == variant) & (per_farm_df['dataset'] == ds)]
        jitter = np.random.default_rng(1).uniform(-0.12, 0.12, len(d))
        ax.scatter(np.full(len(d), i) + jitter, d['bal_acc'], color=DS_COLOURS[ds], s=36, alpha=0.8)
        s = summary[(summary['variant'] == variant) & (summary['dataset'] == ds)].iloc[0]
        ax.hlines(s['pooled_bal_acc'], i - 0.3, i + 0.3, color=DS_COLOURS[ds], lw=2.5)
        ax.text(i, 0.04, f'p {s["perm_p"]:.3f}', ha='center', fontsize=8)
        if ds in old:   # old random-fold score, dashed
            ax.hlines(old[ds], i - 0.3, i + 0.3, color='grey', lw=1.2, ls='--')
    ax.axhline(0.5, color='black', lw=0.8, ls=':')          # chance
    ax.set_xticks(range(len(DATASETS)), DATASETS)
    ax.set_ylim(0, 1)
    ax.set_title(variant, fontsize=11)
axes[0].set_ylabel('Balanced accuracy on the held-out farm')
fig.suptitle('Bagged vs unbagged, leave-one-farm-out: dots = farms, solid = pooled, '
             'dashed = notebook 03 random folds', fontsize=10)
fig.tight_layout()
fig.savefig(FIG_DIR / f'{PREFIX}per_farm{SET_TAG}.png', dpi=FIG_DPI)

print(summary[['dataset', 'variant', 'pooled_bal_acc', 'pooled_roc_auc', 'mean_farm_bal_acc',
               'farms_above_chance', 'perm_p']].round(3).to_string(index=False))
print('03b done.')
