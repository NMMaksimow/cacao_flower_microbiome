"""Build notebooks/03_16S_ITS1_RFC.ipynb - part 4: sensitivity panels and importance."""
import json, pathlib

CELLS = []
def md(src):   CELLS.append(('markdown', src))
def code(src): CELLS.append(('code', src))

# ── Transformation sensitivity ───────────────────────────────────────────────
md("""## 12. Does the transformation drive the result?

CLR is the default here, but the PCoA in 08a/09c uses rCLR, and one could argue for
raw counts. Random Forest is not invariant to the choice: CLR and TSS divide by a
per-sample quantity, which can reorder samples within a feature, and that ordering is
what a split threshold acts on.

The specific risk with raw counts is documented in 08a: unbagged flowers are
systematically sequenced deeper than bagged ones on every farm. On untransformed
counts a forest can reach high accuracy on `sample_type` by reading library size
rather than community composition. The diagnostic below quantifies that directly,
and the panel shows whether it materialises.

Fixed hyperparameters and a single 5-fold CV are used here — this panel compares
transforms against each other, not against an absolute standard.""")

code('''# ── Diagnostic: how strongly is library size tied to each target? ───────────
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
print('\\nHigh eta^2 means raw counts give the forest a technical shortcut.')''')

code('''TRANSFORMS = ['raw', 'tss', 'clr', 'rclr']

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
tr_df.to_csv(TAB_DIR / f'{FIG_PREFIX}transform_comparison.csv', index=False)''')

code('''fig, axes = plt.subplots(2, len(CLF_TARGETS), figsize=(4.1 * len(CLF_TARGETS), 7.4),
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
plt.show()''')

# ── Prevalence sensitivity ───────────────────────────────────────────────────
md("""## 13. Effect of the prevalence threshold

10% is used throughout, to match 08a/08b and the ANCOMBC2 scripts — RF importances are
only comparable with ANCOMBC2 and ISA results if all three see the same feature space.
The cost differs sharply between markers:

| threshold | ITS1 features | ITS1 reads kept | 16S features | 16S reads kept |
|---|---|---|---|---|
| 5 % | 554 | 96.9 % | 580 | 77.6 % |
| **10 %** | **303** | **94.2 %** | **231** | **65.2 %** |
| 20 % | 158 | 87.1 % | 87 | 46.1 % |

For ITS1 the step from 5% to 10% costs under 3% of reads. For 16S it discards a third
of them — the bacterial community is far more even. That number belongs in the Methods.

The last row, `10% within unbagged`, computes prevalence inside the unbagged subset
instead of across all 294 flowers. Three of the targets live only on unbagged flowers,
and a taxon common there but rare overall is dropped by the global filter: 45 ITS1
features (0.6% of reads) and 37 16S features (3.8% of reads).""")

code('''unbagged_ids = list(meta.index[meta['sample_type'] == 'unbagged_flower'])
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
pv_df.to_csv(TAB_DIR / f'{FIG_PREFIX}prevalence_comparison.csv', index=False)''')

code('''fig, axes = plt.subplots(2, len(CLF_TARGETS), figsize=(4.3 * len(CLF_TARGETS), 7.8),
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
plt.show()''')

# ── Fold-averaged importances ────────────────────────────────────────────────
md("""## 14. Permutation and SHAP importance

Both measures are computed **fold-wise on held-out data**, not on a model refitted to
everything. That matters: with `min_samples_leaf=1` a forest fitted on all samples
reaches balanced accuracy 1.000, and shuffling any one of several hundred correlated
features then changes nothing — the first run produced permutation importances of
exactly zero for 8 of 12 combinations. Measuring on the test fold removes that
failure mode.

Error bars are the spread across folds. A tall bar with a wide error bar surfaced in
one fold; a tall bar with a tight one is a taxon the model leans on consistently.""")

code('''def plot_importance_panel(measure, title, xlabel, fname):
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
            err = imp[f'{measure}_std'].reindex(top.index)
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
pd.concat({k: v for k, v in
           ((f'{t}|{d}|{v}', RESULTS[(t, d, v)]['imp']) for (t, d, v) in RESULTS)},
          names=['target_dataset_variant', 'feature']) \\
  .reset_index().to_csv(TAB_DIR / f'{FIG_PREFIX}importances.csv', index=False)
print('Saved importances table.')''')

# ── Taxonomic importance heatmap ─────────────────────────────────────────────
md("""## 15. Importance aggregated by taxonomic rank

Gini importance summed per taxon, using the deepest rank each feature was classified
to: `g__Fusarium` when the genus is known, otherwise `f__Nectriaceae`, `o__Hypocreales`
and so on. That fallback replaces the former single `unclassified` bin, which pooled
55 of 303 ITS1 features and 47 of 231 16S features and therefore always topped the
chart — its height was an artefact of aggregation, not one important taxon.

Prefixes are kept so the labels match the ANCOMBC2 and ISA figures.

**Columns are normalised independently** (each target's maximum becomes 1), so colour
compares taxa *within* a column. A bright cell under a weak target is the maximum of
near-zero values and does not mean the same absolute importance as a bright cell under
a strong one.""")

code('''def imp_by_rank(imp_series, ids, rank_map):
    """Sum fold-averaged Gini importance per deepest-assigned taxonomic rank."""
    s = pd.Series(imp_series.values, index=ids)
    key = pd.Series({i: rank_map.get(i, str(i)[:12]) for i in ids})
    return s.groupby(key).sum().sort_values(ascending=False)


def rank_heatmap(ds, n_top=25):
    cols = {}
    for target in CLF_TARGETS:
        r = RESULTS[(target, ds, 'raw')]
        cols[target] = imp_by_rank(r['imp']['gini_mean'],
                                   FEATURES[ds]['ids'], FEATURES[ds]['rank'])
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
    rank_heatmap(ds).to_csv(TAB_DIR / f'{FIG_PREFIX}rank_importance_{ds}.csv')''')

pathlib.Path(__file__).with_name('cells_part4.json').write_text(
    json.dumps(CELLS, ensure_ascii=False))
print(f'part 4: {len(CELLS)} cells')
