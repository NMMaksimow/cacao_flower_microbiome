#!/usr/bin/env python
# coding: utf-8
"""03b vis - leave-one-farm-out RFC: predicted probability for every flower.

Reads the held-out predictions written by 03b_16S_ITS1_RFC_lofo.py; fits nothing.

For each held-out farm, every flower is a point at its predicted probability of
being unbagged, in two columns: flowers that were truly bagged (left) and truly
unbagged (right). A model that separates the treatments on a farm it never saw
puts the left column low and the right column high.

Batch check: the flowers of extraction batch 2024-09-20 (the only batch holding
both treatments, mt and yb) are drawn as black-edged diamonds. If the model had
learned extraction batches rather than treatment, these flowers, extracted
together, would not sort by treatment.

Run from the notebooks/ folder:
  CFM_SAMPLE_SET=main python 03b_16S_ITS1_RFC_lofo_vis.py
Output: results/<version>/rfc_lofo/figures/03b_lofo_probabilities_<set>.png
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# ── Parameters ───────────────────────────────────────────────────────────────
CFM_VERSION = os.environ.get('CFM_VERSION', 'v2')
SAMPLE_SET  = os.environ.get('CFM_SAMPLE_SET', 'main')
SET_TAG     = '_288' if SAMPLE_SET == 'main' else '_294'
LOFO_DIR    = Path('..').resolve() / 'results' / CFM_VERSION / 'rfc_lofo'

VARIANT     = 'raw'          # 'raw' or 'farm-centred': which LOFO variant to draw
BATCH_CHECK = '2024-09-20'   # extraction batch holding both treatments
FARM_LEVELS = ['ib', 'vr', 'sa', 'kk', 'mt', 'vi', 'yb']   # inside forest -> full sun
DATASETS    = ['ITS1', '16S', 'Combined']
DS_COLOURS  = {'ITS1': '#9467BD', '16S': '#D62728', 'Combined': '#2CA02C'}  # as notebook 03

POINT_SIZE  = 14     # ordinary flowers
POINT_ALPHA = 0.55   # transparency of ordinary flowers (overlap is common)
BATCH_SIZE  = 42     # diamonds of the batch-check flowers
JITTER      = 0.13   # horizontal spread of the points within a column
COL_OFFSET  = 0.2    # distance of the bagged / unbagged columns from the farm centre
FIG_SIZE    = (13, 9.5)   # inches, three rows of panels
FIG_DPI     = 300

# ── Data ─────────────────────────────────────────────────────────────────────
pred = pd.read_csv(LOFO_DIR / 'tables' / f'03b_lofo_predictions{SET_TAG}.csv')
per_farm = pd.read_csv(LOFO_DIR / 'tables' / f'03b_lofo_per_farm{SET_TAG}.csv')
pred = pred[pred['variant'] == VARIANT]
per_farm = per_farm[per_farm['variant'] == VARIANT]

rng = np.random.default_rng(1)   # fixed jitter, so reruns give the same figure

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(len(DATASETS), 1, figsize=FIG_SIZE, sharex=True, sharey=True)
for ax, ds in zip(axes, DATASETS):
    d = pred[pred['dataset'] == ds]
    for i, farm in enumerate(FARM_LEVELS):
        for true, offset in (('bagged_flower', -COL_OFFSET), ('unbagged_flower', COL_OFFSET)):
            f = d[(d['farm_id'] == farm) & (d['true'] == true)]
            x = i + offset + rng.uniform(-JITTER, JITTER, len(f))
            batch = (f['extraction_batch'] == BATCH_CHECK).values
            # Bagged flowers as open circles, unbagged as filled circles
            face = 'none' if true == 'bagged_flower' else DS_COLOURS[ds]
            ax.scatter(x[~batch], f['proba_unbagged'].values[~batch], s=POINT_SIZE,
                       facecolors=face, edgecolors=DS_COLOURS[ds], alpha=POINT_ALPHA, lw=0.8)
            ax.scatter(x[batch], f['proba_unbagged'].values[batch], s=BATCH_SIZE, marker='D',
                       facecolors=face, edgecolors='black', lw=1.0, zorder=3)
        # Balanced accuracy on this held-out farm, above its two columns
        acc = per_farm[(per_farm['dataset'] == ds) & (per_farm['farm'] == farm)]['bal_acc'].iloc[0]
        ax.text(i, 1.04, f'{acc:.2f}', ha='center', va='bottom', fontsize=8, color='grey')
    ax.axhline(0.5, color='black', lw=0.8, ls=':')    # decision threshold
    for i in range(1, len(FARM_LEVELS)):              # separators between farms
        ax.axvline(i - 0.5, color='grey', lw=0.4, alpha=0.5)
    ax.set_ylim(-0.03, 1.13)
    ax.set_ylabel(f'{ds}\nP(unbagged)')
    ax.yaxis.label.set_color(DS_COLOURS[ds])

axes[-1].set_xticks(range(len(FARM_LEVELS)), FARM_LEVELS)
axes[-1].set_xlabel('Held-out farm (left column: truly bagged, right column: truly unbagged; '
                    'grey number: balanced accuracy on that farm)')

legend = [
    Line2D([], [], marker='o', ls='', mfc='none', mec='grey', label='Bagged flower'),
    Line2D([], [], marker='o', ls='', mfc='grey', mec='grey', label='Unbagged flower'),
    Line2D([], [], marker='D', ls='', mfc='none', mec='black', label=f'Extraction batch {BATCH_CHECK}'),
    Line2D([], [], color='black', ls=':', lw=0.8, label='Decision threshold 0.5'),
]
axes[0].legend(handles=legend, loc='lower center', bbox_to_anchor=(0.5, 1.12), ncol=4,
               frameon=False, fontsize=9)   # legend above the top panel
fig.suptitle(f'Leave-one-farm-out RFC, {VARIANT} features: predicted probability of being '
             f'unbagged for every flower ({SET_TAG[1:]} set)', fontsize=11, y=0.995)
fig.tight_layout()
out = LOFO_DIR / 'figures' / f'03b_lofo_probabilities{SET_TAG}.png'
fig.savefig(out, dpi=FIG_DPI, facecolor='white')
print(f'Saved {out}')
