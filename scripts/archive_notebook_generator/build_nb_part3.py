"""Build notebooks/03_16S_ITS1_RFC.ipynb - part 3: summary, figures, sensitivity panels."""
import json, pathlib

CELLS = []
def md(src):   CELLS.append(('markdown', src))
def code(src): CELLS.append(('code', src))

# ── Summary ──────────────────────────────────────────────────────────────────
md("""## 7. Summary of all models

`nested_bal_acc` is the headline, unbiased estimate. `perm_observed` comes from the
permutation procedure (fixed hyperparameters, single CV) and is only meaningful next
to its own p-value — the two columns answer different questions and are not expected
to be identical.""")

code('''rows = []
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
    show[c] = show[c].map(lambda v: '\u2014' if pd.isna(v) else f'{v:.3f}')
show['perm_p'] = show['perm_p'].map(lambda v: f'{v:.4f}')
print(show[['target', 'dataset', 'variant', 'n', 'chance', 'nested_bal_acc',
            'nested_std', 'above_chance', 'nested_mcc', 'roc_auc',
            'perm_p', 'n_perm']].to_string(index=False))

print('\\nRegression (pi_clsm):')
for (ds, variant), r in REG_RESULTS.items():
    print(f'  {ds:9s} {variant:12s} n={r["n"]}  Spearman={r["spearman"]:.3f}  '
          f'R2={r["r2"]:.3f}  MAE={r["mae"]:.1f}  p={r["perm"]["pvalue"]:.4f}')''')

code('''# ── Heatmap: how far each model gets from chance towards perfect ────────────
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
ax.set_title('Nested-CV performance above chance\\n'
             'cell annotation = raw balanced accuracy', fontweight='bold', fontsize=11)
ax.set_xlabel(''); ax.set_ylabel('')
plt.tight_layout()
fig.savefig(FIG_DIR / f'{FIG_PREFIX}summary_heatmap.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()''')

# ── Confusion matrices + Gini ────────────────────────────────────────────────
md("""## 8. Confusion matrices and Gini importance

One figure per target, one column per feature set. Confusion matrices come from the
out-of-fold predictions (fixed hyperparameters), so no sample is scored by a model
that saw it during training.""")

code('''def plot_target_panel(target, variant='raw', n_top=N_TOP_IMP):
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
            f'{ds}  (n={r["n"]})\\n'
            f'nested bal.acc = {r["nested"]["bal_acc_mean"]:.3f} '
            f'± {r["nested"]["bal_acc_std"]:.3f}\\n'
            f'MCC = {r["nested"]["mcc_mean"]:.3f}   ROC-AUC = {roc}\\n'
            f'permutation p = {r["perm"]["pvalue"]:.4f} '
            f'(n={r["perm"].get("n_perm", N_PERM)})', fontsize=9)
        axes[0, col].tick_params(axis='x', rotation=30, labelsize=8)

        imp = r['imp']['gini_mean'].nlargest(n_top).sort_values()
        err = r['imp']['gini_std'].reindex(imp.index)
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
        plot_target_panel(t, v)''')

# ── Permutation tests ────────────────────────────────────────────────────────
md("""## 9. Permutation tests

Each histogram is the distribution of balanced accuracy obtained after shuffling the
labels, i.e. what this design, this sample size and this class imbalance produce when
there is no signal. The dashed line is the observed score. A fixed chance line of
1/n_classes is the theoretical expectation; the empirical null shows the spread around
it, which is what a p-value needs.""")

code('''perm_rows = [(t, v) for t in CLF_TARGETS
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
        ax.set_title(f'{target}{tag} | {ds}\\n'
                     f'observed = {r["perm"]["observed"]:.3f}, '
                     f'p = {r["perm"]["pvalue"]:.4f} '
                     f'({r["perm"].get("n_perm", N_PERM)} shuffles)', fontsize=8.5)
        ax.set_xlabel('Balanced accuracy'); ax.set_ylabel('Permutations')
        ax.tick_params(labelsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(FIG_DIR / f'{FIG_PREFIX}permutation_tests.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()''')

# ── Calibration ──────────────────────────────────────────────────────────────
md("""## 10. Probability calibration — binary targets

Random Forest probabilities on small samples tend to be pushed toward 0 and 1. This
matters whenever `predict_proba` is read as a probability rather than a ranking. The
diagonal is perfect calibration; the Brier score summarises the gap (lower is better).""")

code('''BINARY_TARGETS = [t for t in CLF_TARGETS if len(TARGETS[t]['order'] or []) == 2]

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
    ax.set_title(f'{target}\\npositive class = '
                 f'"{RESULTS[(target, DATASETS[0], "raw")]["classes"][1]}"', fontsize=10)
    ax.set_xlabel('Mean predicted probability'); ax.set_ylabel('Observed frequency')
    ax.legend(fontsize=8, loc='upper left')

plt.tight_layout()
fig.savefig(FIG_DIR / f'{FIG_PREFIX}calibration.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()''')

# ── Regression figure ────────────────────────────────────────────────────────
md("""## 11. Regression diagnostics — `pi_clsm`

Out-of-fold predictions against observed pollen load. Points are coloured by the
3-class grouping used in section 5, so the two framings can be read against each other.""")

code('''POLL_COLOURS = {'none': '#4C72B0', 'weak': '#DD8452', 'strong': '#55A868'}

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

        ax.set_title(f'{ds} | {variant}  (n={r["n"]})\\n'
                     f'Spearman rho = {r["spearman"]:.3f}  '
                     f'(p = {r["perm"]["pvalue"]:.4f})\\n'
                     f'R² = {r["r2"]:.3f}   MAE = {r["mae"]:.1f}', fontsize=9)
        ax.set_xlabel(xlab)
        if col == 0:
            ax.set_ylabel('Predicted (out-of-fold)')
        if row == 0 and col == 0:
            ax.legend(title='3-class group', fontsize=8, title_fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(FIG_DIR / f'{FIG_PREFIX}regression_pi_clsm.png', dpi=FIG_DPI,
            bbox_inches='tight', facecolor='white')
plt.show()''')

pathlib.Path(__file__).with_name('cells_part3.json').write_text(
    json.dumps(CELLS, ensure_ascii=False))
print(f'part 3: {len(CELLS)} cells')
