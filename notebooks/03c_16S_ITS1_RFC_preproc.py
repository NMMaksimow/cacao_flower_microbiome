#!/usr/bin/env python
# coding: utf-8
"""03c - Random Forest: does the preprocessing drive the result?

Two questions that one 5-fold CV run (former sections 12-13 of notebook 03) could not
settle, because its scores move by several points from one fold split to the next:

  A. Feature representation. Why do raw counts and TSS score higher than CLR?
  B. Prevalence threshold. Does the 10% filter remove taxa of small classes (strongly
     pollinated flowers: 15; fungal score 5: 20)?

Design
  Model     as score_fixed() of notebook 03: 500 trees, max_features='sqrt',
            class_weight='balanced', no tuning, no farm centring
  CV        RepeatedStratifiedKFold 5 folds x 5 repeats, the same splits for every
            setting (paired comparison); one score per repeat = mean of its 5 folds
  Targets   farm_id, sample_type (all flowers); pollination, fungal_colonization
            (unbagged flowers), defined as in notebook 03
  Markers   ITS1 OTU97, 16S ASV
  Block A   representations, all on the 10%-of-all-flowers feature set:
              raw         counts
              tss         counts / sample total
              clr_czm     CLR from R (08c), CZM zero replacement -- used in notebook 03
              clr_delta   CLR with the fixed delta = 1/D^2 used before 2026-11
                          (kept here only to show what the switch changed)
              rclr        log of non-zero values minus their sample mean; zeros stay 0
              presence    1 if the taxon has reads, else 0
              log_depth   one feature: log of the sample's total reads (control: how
                          far does sequencing depth alone get?)
  Block B   prevalence threshold, TSS and CLR-CZM:
              all00 (TSS only), all05, all10, all20  >= x% of all flowers
              unbagged05, unbagged10                 >= x% of unbagged flowers
            08c has no 0% CLR matrix (cmultRepl on ~26 700 16S features takes ~45 min);
            TSS needs no zero replacement, so 0% is tested on TSS.
  Null      pollination and fungal_colonization, clr_czm all10: labels shuffled
            100 times, each scored with the same 5 x 5 procedure. The null shows how
            high a score this procedure reaches with no signal at all.
  Diagnostics (cheap)
            eta^2 of log sequencing depth for each target; for each representation,
            the share of (zero, non-zero) pairs within a feature where the zero ranks
            ABOVE the non-zero. 0 means 'taxon absent' stays a clean split for a tree.

Every finished score is appended to the CSV at once; a restarted run skips what is
already there. Delete the CSVs to recompute from scratch.

Environment
  CFM_VERSION        v1 | v2 (default v2)
  CFM_SAMPLE_SET     main (288, default) | with_deviations (294; needs the 08c
                     matrices of that set)
  RFC_TAG            tag in the output names (03c_hpc_ on Hamilton), as notebook 03
  RFC_PREPROC_QUICK  1 = quick pass to test the code (100 trees, 2 repeats, 5 shuffles)

Run from the notebooks/ folder, on Hamilton via scripts/03c_hpc_rfc_preproc.sh.
Figures: 03c_16S_ITS1_RFC_preproc_vis.py (reads the tables, fits nothing).
"""

import os
import time
import warnings
from pathlib import Path

import biom
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score

from utils_rfc import SET_TAG, load_clr   # CLR matrices of scripts/08c_16S_ITS1_clr_export.R

warnings.filterwarnings('ignore')

# ── 0. Parameters ────────────────────────────────────────────────────────────
CFM_VERSION = os.environ.get('CFM_VERSION', 'v2')
SAMPLE_SET  = os.environ.get('CFM_SAMPLE_SET', 'main')
OUT_TAG     = os.environ.get('RFC_TAG', '')
QUICK       = os.environ.get('RFC_PREPROC_QUICK', '0') == '1'
assert CFM_VERSION in ('v1', 'v2') and SAMPLE_SET in ('main', 'with_deviations')

ROOT        = Path('..').resolve()
VER_DIR     = ROOT / 'results' / CFM_VERSION
BIOM_DIR    = VER_DIR / 'biom'
BIOM_SUFFIX = '' if SAMPLE_SET == 'main' else '_with_deviations'
TAB_DIR     = VER_DIR / 'tables'
TAB_DIR.mkdir(parents=True, exist_ok=True)
PREFIX      = f'03c{OUT_TAG}_' + ('quick_' if QUICK else '')
SCORES_CSV  = TAB_DIR / f'{PREFIX}scores{SET_TAG[SAMPLE_SET]}.csv'
NULL_CSV    = TAB_DIR / f'{PREFIX}null{SET_TAG[SAMPLE_SET]}.csv'
DIAG_CSV    = TAB_DIR / f'{PREFIX}diagnostics{SET_TAG[SAMPLE_SET]}.csv'

RANDOM_STATE   = 42
N_SPLITS       = 5
PI_THRESH_POLL = 40     # pollen grains per stigma; none / weak / strong, as notebook 03
if QUICK:
    N_TREES, N_REPEATS, N_NULL = 100, 2, 5
else:
    N_TREES, N_REPEATS, N_NULL = 500, 5, 100

MARKERS = {'ITS1': 'cfm_its1_otu97_fungi_biosamples', '16S': 'cfm_16s_bacteria_biosamples'}
NULL_TARGETS = ['pollination', 'fungal_colonization']

print(f'Profile {"QUICK" if QUICK else "FULL"} | data {CFM_VERSION}, set {SAMPLE_SET} | '
      f'trees {N_TREES}, CV {N_SPLITS} x {N_REPEATS}, null {N_NULL} -> {SCORES_CSV.name}')


# ── 1. Data ──────────────────────────────────────────────────────────────────
def load_counts(stem):
    """features x samples count table from the biom export of 01a/01b."""
    return biom.load_table(str(BIOM_DIR / f'{stem}{BIOM_SUFFIX}_feature-table.biom')) \
               .to_dataframe(dense=True)


meta = pd.read_csv(ROOT / 'data' / 'sample_metadata_with_microscopy.txt',
                   sep='\t', index_col=0, dtype=str)
meta = meta.drop(index='#q2:types', errors='ignore')
meta = meta[meta['sample_type'].isin(['bagged_flower', 'unbagged_flower'])]
for col in ['pi_clsm', 'fungal_colonization_score']:
    meta[col] = pd.to_numeric(meta[col].astype(str).str.replace(',', '.', regex=False),
                              errors='coerce')


def make_targets(m):
    """The four classification targets of notebook 03, as label Series."""
    unb = m[m['sample_type'] == 'unbagged_flower']
    pol = unb['pi_clsm'].dropna().apply(
        lambda v: 'none' if v == 0 else ('weak' if v < PI_THRESH_POLL else 'strong'))
    fc = unb['fungal_colonization_score'].dropna()
    fc = fc[fc.isin([0, 5])].map({0.0: 'none', 5.0: 'strong'})
    return {'farm_id': m['farm_id'], 'sample_type': m['sample_type'],
            'pollination': pol, 'fungal_colonization': fc}


# ── 2. Representations ───────────────────────────────────────────────────────
# Input `c`: counts, samples x features. Each function returns samples x features.
def tss(c):
    return c.div(c.sum(axis=1), axis=0)


def clr_delta(c):
    """CLR with the fixed replacement delta = 1/D^2 used before 2026-11 (comparison only)."""
    ra = tss(c).values
    delta = 1.0 / ra.shape[1] ** 2
    k = (ra == 0).sum(axis=1, keepdims=True)          # zeros per sample
    comp = np.where(ra == 0, delta, ra * (1 - k * delta))
    logc = np.log(comp)
    return pd.DataFrame(logc - logc.mean(axis=1, keepdims=True), index=c.index, columns=c.columns)


def rclr(c):
    """Robust CLR as vegan: mean log over the non-zero values; zeros stay 0."""
    ra = tss(c).values
    nz = ra > 0
    logx = np.where(nz, np.log(np.where(nz, ra, 1.0)), 0.0)
    gm = logx.sum(axis=1) / np.maximum(nz.sum(axis=1), 1)
    return pd.DataFrame(np.where(nz, logx - gm[:, None], 0.0), index=c.index, columns=c.columns)


def zero_rank_overlap(X, present):
    """Mean over features of P(a zero ranks above a non-zero of the same feature).

    0 = presence / absence is a clean split for a tree; larger = absent and present
    samples are interleaved. X, present: samples x features, same order.
    """
    X, P = X.values, present.values
    shares = []
    for j in range(X.shape[1]):
        z, nz = X[~P[:, j], j], X[P[:, j], j]
        if len(z) and len(nz):
            shares.append((z[:, None] > nz[None, :]).mean())
    return float(np.mean(shares))


def eta2(values, groups):
    """Share of the variance of `values` explained by `groups` (one-way layout)."""
    grand = values.mean()
    between = sum(len(v) * (v.mean() - grand) ** 2 for _, v in values.groupby(groups))
    return float(between / ((values - grand) ** 2).sum())


# ── 3. Scoring with incremental saving ───────────────────────────────────────
def make_rf():
    # n_jobs=1: parallelism lives in the CV loop (as notebook 03)
    return RandomForestClassifier(n_estimators=N_TREES, max_features='sqrt',
                                  class_weight='balanced', random_state=RANDOM_STATE,
                                  n_jobs=1)


def repeat_scores(X, y):
    """Balanced accuracy, one value per repeat (mean over its 5 folds)."""
    cv = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS,
                                 random_state=RANDOM_STATE)
    s = cross_val_score(make_rf(), X.loc[y.index].values, y.values, cv=cv,
                        scoring='balanced_accuracy', n_jobs=-1)
    return s.reshape(N_REPEATS, N_SPLITS).mean(axis=1)


def done_keys(path, cols):
    """Keys already in a results CSV, so a restarted run can skip them."""
    if not path.exists():
        return set()
    return set(map(tuple, pd.read_csv(path, dtype=str)[cols].values))


def append_row(path, row):
    pd.DataFrame([row]).to_csv(path, mode='a', header=not path.exists(), index=False)


SCORE_KEY = ['block', 'marker', 'features', 'config', 'target']
done = done_keys(SCORES_CSV, SCORE_KEY)
t0 = time.time()


def score(block, marker, features, config, target, X, y):
    key = (block, marker, features, config, target)
    if key in done:
        return
    y = y.loc[y.index.intersection(X.index)]
    s = repeat_scores(X, y)
    append_row(SCORES_CSV, dict(zip(SCORE_KEY, key), n=len(y), n_features=X.shape[1],
                                mean=s.mean(), sd=s.std(ddof=1) if len(s) > 1 else np.nan,
                                min=s.min(), max=s.max(),
                                repeats=';'.join(f'{v:.4f}' for v in s)))
    print(f'  {block} {marker:4s} {features:9s} {config:10s} {target:19s} '
          f'{s.mean():.3f} +- {s.std():.3f}  [{time.time() - t0:6.0f} s]', flush=True)


# ── 4. Run ───────────────────────────────────────────────────────────────────
diag_rows = []
for marker, stem in MARKERS.items():
    counts_all = load_counts(stem)                              # features x samples
    clr = {cfg: load_clr(marker, VER_DIR, SAMPLE_SET, cfg, counts=counts_all)
           for cfg in ['all10', 'all05', 'all20', 'unbagged05', 'unbagged10']}
    samples = clr['all10'].index.intersection(meta.index)
    counts = counts_all.T.loc[samples]                          # samples x features
    targets = make_targets(meta.loc[samples])
    print(f'{marker}: {len(samples)} flowers, {counts.shape[1]} features in the biom')

    # Feature set of each prevalence config = the columns of its 08c matrix, so TSS
    # and CLR always see exactly the same taxa. all00 = every taxon with reads.
    feats = {cfg: X.columns for cfg, X in clr.items()}
    feats['all00'] = counts.columns[counts.sum(axis=0) > 0]

    c10 = counts[feats['all10']]
    reps = {
        'raw':       c10,
        'tss':       tss(c10),
        'clr_czm':   clr['all10'].loc[samples],
        'clr_delta': clr_delta(c10),
        'rclr':      rclr(c10),
        'presence':  (c10 > 0).astype(int),
        'log_depth': np.log(counts.sum(axis=1)).to_frame('log_depth'),
    }

    # Diagnostics (seconds)
    for t, y in targets.items():
        d = np.log(counts.loc[y.index].sum(axis=1))
        diag_rows.append(dict(marker=marker, measure='eta2_log_depth', item=t,
                              value=eta2(d, y.values)))
    for name in ['tss', 'clr_czm', 'clr_delta', 'rclr']:
        diag_rows.append(dict(marker=marker, measure='zero_rank_overlap', item=name,
                              value=zero_rank_overlap(reps[name], c10 > 0)))
    pd.DataFrame(diag_rows).to_csv(DIAG_CSV, index=False)

    # Block A: representations
    for name, X in reps.items():
        for t, y in targets.items():
            score('A', marker, name, 'all10', t, X, y)

    # Block B: prevalence thresholds (all10 is shared with block A but scored again
    # here, so each block is complete on its own)
    for cfg in ['all00', 'all05', 'all10', 'all20', 'unbagged05', 'unbagged10']:
        for t, y in targets.items():
            score('B', marker, 'tss', cfg, t, tss(counts[feats[cfg]]), y)
            if cfg in clr:
                score('B', marker, 'clr_czm', cfg, t, clr[cfg].loc[samples], y)

    # Null: shuffled labels, same procedure, clr_czm all10
    null_done = done_keys(NULL_CSV, ['marker', 'target', 'perm'])
    for t in NULL_TARGETS:
        y = targets[t]
        for i in range(N_NULL):
            if (marker, t, str(i)) in null_done:
                continue
            y_perm = pd.Series(np.random.default_rng(i).permutation(y.values), index=y.index)
            s = repeat_scores(reps['clr_czm'], y_perm)
            append_row(NULL_CSV, dict(marker=marker, target=t, perm=i, mean=s.mean()))
        print(f'  null {marker} {t}: {N_NULL} shuffles  [{time.time() - t0:6.0f} s]', flush=True)

print(f'done in {time.time() - t0:.0f} s')
