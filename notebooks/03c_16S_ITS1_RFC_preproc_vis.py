#!/usr/bin/env python
# coding: utf-8
"""03c vis - figures of the preprocessing check; reads the tables, fits nothing.

  03c<tag>_transforms<set>.png   block A: balanced accuracy per feature representation
  03c<tag>_prevalence<set>.png   block B: balanced accuracy per prevalence threshold

Bars: mean over the 5 repeats of 5-fold CV; whiskers: SD over repeats; dots: the
individual repeats. Dotted line: chance (1 / number of classes). Dashed red line
(pollination, fungal_colonization): 95th percentile of the shuffled-label null, i.e.
the score this procedure reaches without any signal in 1 run of 20.

Run from the notebooks/ folder, after copying the 03c tables from Hamilton:
  python 03c_16S_ITS1_RFC_preproc_vis.py
Same environment variables as 03c_16S_ITS1_RFC_preproc.py (CFM_VERSION,
CFM_SAMPLE_SET, RFC_TAG, RFC_PREPROC_QUICK); RFC_TAG defaults to _hpc here,
because the tables are produced on Hamilton.
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils_rfc import SET_TAG

# ── Parameters ───────────────────────────────────────────────────────────────
CFM_VERSION = os.environ.get('CFM_VERSION', 'v2')
SAMPLE_SET  = os.environ.get('CFM_SAMPLE_SET', 'main')
OUT_TAG     = os.environ.get('RFC_TAG', '_hpc')
QUICK       = os.environ.get('RFC_PREPROC_QUICK', '0') == '1'
VER_DIR     = Path('..').resolve() / 'results' / CFM_VERSION
TAB_DIR     = VER_DIR / 'tables'
FIG_DIR     = VER_DIR / 'figures' / 'random_forest'
FIG_DIR.mkdir(parents=True, exist_ok=True)
PREFIX      = f'03c{OUT_TAG}_' + ('quick_' if QUICK else '')
TAG         = SET_TAG[SAMPLE_SET]

MARKERS  = ['ITS1', '16S']                                   # figure rows
TARGETS  = ['farm_id', 'sample_type', 'pollination', 'fungal_colonization']   # columns
N_CLASSES = {'farm_id': 7, 'sample_type': 2, 'pollination': 3, 'fungal_colonization': 2}

# Block A: bar order and labels on the x axis
REPS = ['raw', 'tss', 'clr_czm', 'clr_delta', 'rclr', 'presence', 'log_depth']
REP_LABELS = {'raw': 'raw', 'tss': 'TSS', 'clr_czm': 'CLR\nCZM', 'clr_delta': 'CLR\ndelta',
              'rclr': 'rCLR', 'presence': 'pres./\nabs.', 'log_depth': 'depth\nonly'}
# Block B: threshold order and labels
CONFIGS = ['all00', 'all05', 'all10', 'all20', 'unbagged05', 'unbagged10']
CONFIG_LABELS = {'all00': '0%', 'all05': '5%', 'all10': '10%\n(used)', 'all20': '20%',
                 'unbagged05': '5%\nunbagged', 'unbagged10': '10%\nunbagged'}

# Colours: the representation used in notebook 03 (CLR-CZM) in the marker colour of
# notebook 03, every alternative in grey; the depth-only control lighter still.
DS_COLOURS  = {'ITS1': '#9467BD', '16S': '#D62728'}
ALT_COLOUR  = '#BBBBBB'
CTRL_COLOUR = '#E5E5E5'
NULL_COLOUR = '#B2182B'      # dashed line: 95th percentile of the null
FEAT_COLOUR = '#666666'      # block B: dashed line, number of features (right axis)

BAR_WIDTH   = 0.38           # block B: width of each of the two grouped bars
DOT_SIZE    = 9              # individual repeats
PANEL_W, PANEL_H = 4.2, 3.6  # inches per panel
FIG_DPI     = 300

# ── Data ─────────────────────────────────────────────────────────────────────
scores = pd.read_csv(TAB_DIR / f'{PREFIX}scores{TAG}.csv')
null = pd.read_csv(TAB_DIR / f'{PREFIX}null{TAG}.csv')
null95 = null.groupby(['marker', 'target'])['mean'].quantile(0.95)
diag = pd.read_csv(TAB_DIR / f'{PREFIX}diagnostics{TAG}.csv')
print(diag.pivot_table(index=['measure', 'item'], columns='marker', values='value').round(4))


def repeats(row):
    return [float(v) for v in str(row['repeats']).split(';')]


def reference_lines(ax, marker, target):
    """Chance level and, where a null exists, its 95th percentile."""
    ax.axhline(1 / N_CLASSES[target], color='black', ls=':', lw=1)
    if (marker, target) in null95.index:
        ax.axhline(null95[(marker, target)], color=NULL_COLOUR, ls='--', lw=1)


def panel_grid(title):
    fig, axes = plt.subplots(len(MARKERS), len(TARGETS),
                             figsize=(PANEL_W * len(TARGETS), PANEL_H * len(MARKERS)),
                             sharey=True)
    fig.suptitle(title, fontsize=12, fontweight='bold')
    return fig, axes


# ── Figure A: representations ────────────────────────────────────────────────
fig, axes = panel_grid('Feature representation (RF, 5-fold CV x 5 repeats, '
                       '10% prevalence on all flowers)')
a = scores[scores['block'] == 'A']
for r, marker in enumerate(MARKERS):
    for c, target in enumerate(TARGETS):
        ax = axes[r, c]
        sub = a[(a['marker'] == marker) & (a['target'] == target)].set_index('features')
        for i, rep in enumerate(REPS):
            if rep not in sub.index:
                continue
            row = sub.loc[rep]
            colour = (DS_COLOURS[marker] if rep == 'clr_czm'
                      else CTRL_COLOUR if rep == 'log_depth' else ALT_COLOUR)
            ax.bar(i, row['mean'], color=colour, yerr=row['sd'], capsize=2,
                   error_kw=dict(lw=0.8))
            ax.scatter(np.full(len(repeats(row)), i), repeats(row), s=DOT_SIZE,
                       color='black', alpha=0.6, zorder=3)
        reference_lines(ax, marker, target)
        ax.set_xticks(range(len(REPS)))
        ax.set_xticklabels([REP_LABELS[x] for x in REPS], fontsize=7.5)
        ax.set_title(f'{marker} | {target}', fontsize=9, fontweight='bold')
        ax.set_ylim(0, 1.02)
        if c == 0:
            ax.set_ylabel('Balanced accuracy')
plt.tight_layout()
out = FIG_DIR / f'{PREFIX}transforms{TAG}.png'
fig.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white')
print('->', out)

# ── Figure B: prevalence thresholds ──────────────────────────────────────────
fig, axes = panel_grid('Prevalence threshold (RF, 5-fold CV x 5 repeats); '
                       'grey = TSS, colour = CLR-CZM, dashed line = number of features')
b = scores[scores['block'] == 'B']
for r, marker in enumerate(MARKERS):
    for c, target in enumerate(TARGETS):
        ax = axes[r, c]
        sub = b[(b['marker'] == marker) & (b['target'] == target)]
        for k, (feat, colour) in enumerate([('tss', ALT_COLOUR), ('clr_czm', DS_COLOURS[marker])]):
            s = sub[sub['features'] == feat].set_index('config')
            for i, cfg in enumerate(CONFIGS):
                if cfg not in s.index:          # no 0% CLR matrix (see 08c)
                    continue
                x = i + (k - 0.5) * BAR_WIDTH
                ax.bar(x, s.loc[cfg, 'mean'], width=BAR_WIDTH, color=colour,
                       yerr=s.loc[cfg, 'sd'], capsize=2, error_kw=dict(lw=0.8))
                ax.scatter(np.full(len(repeats(s.loc[cfg])), x), repeats(s.loc[cfg]),
                           s=DOT_SIZE, color='black', alpha=0.6, zorder=3)
        reference_lines(ax, marker, target)
        ax.set_xticks(range(len(CONFIGS)))
        ax.set_xticklabels([CONFIG_LABELS[x] for x in CONFIGS], fontsize=7.5)
        ax.set_title(f'{marker} | {target}', fontsize=9, fontweight='bold')
        ax.set_ylim(0, 1.02)
        if c == 0:
            ax.set_ylabel('Balanced accuracy')
        # Number of features on a log axis (0% for 16S is ~26 700, 20% under 100)
        nf = sub[sub['features'] == 'tss'].set_index('config')['n_features'].reindex(CONFIGS)
        ax2 = ax.twinx()
        # no markers on this line, so it cannot be mistaken for the repeat dots
        ax2.plot(range(len(CONFIGS)), nf.values, color=FEAT_COLOUR, ls='--', lw=1)
        ax2.set_yscale('log')
        ax2.tick_params(labelsize=7)
        if c == len(TARGETS) - 1:
            ax2.set_ylabel('# features (log)', fontsize=8)
plt.tight_layout()
out = FIG_DIR / f'{PREFIX}prevalence{TAG}.png'
fig.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor='white')
print('->', out)
