"""Build notebooks/03_16S_ITS1_RFC.ipynb - part 2: targets, nested CV, main loops."""
import json, pathlib

CELLS = []
def md(src):   CELLS.append(('markdown', src))
def code(src): CELLS.append(('code', src))

# ── Targets ──────────────────────────────────────────────────────────────────
md("""## 3. Target definitions

Every target is a function `metadata -> y`, so the sample set is derived from the
metadata rather than hard-coded. The feature matrix is intersected with `y.index`,
which is what makes the unbagged-only targets drop to ~100 samples automatically.""")

code('''def target_farm_id(m):
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
print(counts_df.to_string(index=False))''')

# ── Nested CV machinery ──────────────────────────────────────────────────────
md("""## 4. Cross-validation machinery

**Two numbers per model, and they must not be confused.**

*Nested CV score* — the outer loop measures performance, the inner `GridSearchCV`
picks hyperparameters inside each outer fold. Because selection never sees the outer
test fold, this is the unbiased estimate of generalisation. **This is the number that
goes into the manuscript.**

*Permutation-test observed score* — `permutation_test_score` runs a single CV with one
fixed hyperparameter set and compares it against the same procedure applied to shuffled
labels. The p-value is valid precisely because observed and null go through an identical
procedure — but that procedure is not the nested one, so the observed score here will
differ slightly from the nested-CV score. That is expected, not a discrepancy.

The out-of-fold predictions used for confusion matrices and calibration also use the
fixed hyperparameters: a repeated CV predicts each sample several times, whereas a
confusion matrix needs exactly one prediction per sample.""")

code('''def make_clf(**overrides):
    params = dict(n_estimators=N_TREES, max_features='sqrt', class_weight='balanced',
                  random_state=RANDOM_STATE, n_jobs=-1)
    params.update(overrides)
    return RandomForestClassifier(**params)


def make_reg(**overrides):
    params = dict(n_estimators=N_TREES, max_features='sqrt',
                  random_state=RANDOM_STATE, n_jobs=-1)
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

def cache_load(key):
    p = _cache_path(key)
    if p.exists() and not FORCE_RECOMPUTE:
        return json.loads(p.read_text())
    return None

def cache_save(key, obj):
    _cache_path(key).write_text(json.dumps(obj, indent=1))


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
    cached = cache_load(key)
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
    cache_save(key, out)
    print(f'  {key:42s} {out["runtime_s"]:6.1f}s  params={out["best_params"]}')
    return out


def oof_predictions(X, y, best_params, blocked=False):
    """Single StratifiedKFold pass with fixed hyperparameters: every sample predicted
    exactly once, as a confusion matrix and a calibration curve require. A repeated CV
    predicts each sample several times and cannot serve that purpose."""
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    est, _ = _estimator_and_grid('clf', blocked, **best_params)
    y_pred  = cross_val_predict(est, X, y, cv=cv, method='predict')
    y_proba = cross_val_predict(est, X, y, cv=cv, method='predict_proba')
    return y_pred, y_proba


def run_permutation_test(X, y, best_params, key, kind='clf', blocked=False,
                         n_perm=None):
    """Empirical null from shuffled labels -> p-value for the observed score.

    permutation_test_score puts the observed score and every null replicate through
    the SAME procedure, which is what makes the p-value valid. Fewer trees are used
    here than in the main models (N_TREES_PERM) -- allowed for the same reason.
    """
    cached = cache_load(key)
    if cached is not None:
        print(f'  [cache] {key}')
        return cached

    n_perm = n_perm or N_PERM
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
    cache_save(key, out)
    print(f'  {key:42s} {out["runtime_s"]:6.1f}s  observed={observed:.3f}  '
          f'p={pval:.4f} (n={n_perm})')
    return out


def fold_importances(X, y, best_params, feat_names, blocked=False):
    """All three importance measures, computed fold-wise instead of on a refit model.

    The previous version fitted one model on all data and measured importance there.
    With min_samples_leaf=1 that model memorises the training set (balanced accuracy
    1.000), so shuffling any single feature changes nothing -- hundreds of correlated
    features offer redundant paths to the same leaf, and permutation importance came
    out exactly zero. Measuring on held-out folds removes that failure mode.

    Returns a DataFrame indexed by feature with mean and std across folds for each
    measure; the std separates taxa that matter consistently from one-fold flukes.
    """
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
    return pd.DataFrame(out, index=feat_names)


def spearman_scorer(estimator, X, y):
    """Rank correlation between prediction and truth. Primary regression metric:
    pi_clsm is zero-inflated and right-skewed, so a rank measure is more honest
    than R^2 here."""
    return spearmanr(y, estimator.predict(X))[0]''')

# ── Main classification loop ─────────────────────────────────────────────────
md("""## 5. Fit all classification targets

4 targets x 3 feature sets. Each combination gets nested CV, out-of-fold predictions,
a permutation test, and fold-wise importances.

**Three targets are additionally fitted `within-farm`.** The farm signal dominates this
dataset — ITS1 identifies the farm with ROC-AUC 0.965 — and farm is associated with both
weak biological targets (Cramér V = 0.485 for `fungal_colonization`, 0.326 for
`pollination`). A model could therefore "predict" colonisation by recognising the farm.
Subtracting each farm's mean from every feature leaves only within-farm variation, so
the pair of results answers how much of the signal was farm-mediated.

`farm_id` is not blocked, for obvious reasons.""")

code('''# Which targets get a second, farm-blocked fit.
BLOCKED_TARGETS = ['sample_type', 'pollination', 'fungal_colonization']
VARIANTS = ['raw', 'within-farm']

RESULTS = {}

for target in CLF_TARGETS:
    print(f'\\n=== {target} ===')
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

            n_perm = N_PERM_STRONG if (target, ds) in STRONG_COMBOS and not blocked \\
                     else N_PERM

            nested = nested_cv(X, y, key + '__nested', kind='clf', blocked=blocked)
            bp     = nested['best_params']
            y_pred, y_proba = oof_predictions(X, y, bp, blocked=blocked)
            perm   = run_permutation_test(X, y, bp, key + '__perm', kind='clf',
                                          blocked=blocked, n_perm=n_perm)
            imp    = fold_importances(X, y, bp, FEATURES[ds]['feats'], blocked=blocked)

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
print('\\nDone.')''')

# ── Regression ───────────────────────────────────────────────────────────────
md("""## 6. Target 5 — `pi_clsm` as regression

Thresholding a continuous measurement into classes throws information away and makes
the result depend on where the threshold sits. Regressing `pi_clsm` directly avoids
both, and uses all 102 measured unbagged flowers.

Primary metric is **Spearman rho** between out-of-fold predictions and observations:
`pi_clsm` is zero-inflated (38 of 102 are exactly 0) and right-skewed, so a rank
correlation describes the fit more honestly than R². R² and MAE are reported alongside.

This is not an independent confirmation of the 3-class model in section 5 — it is the
same data and the same question in a different form. Read the two together.

The `within-farm` variant centers **both** the features and the response by farm, which
is the classic within (fixed-effects) estimator. The question becomes: inside one farm,
do flowers richer in a given taxon carry more pollen? Since the raw regression already
sits at rho ~ 0, there is nothing to lose — and anything that appears after centering
would be a far more convincing result than the current one.""")

code('''REG_RESULTS = {}
print('=== pi_clsm (regression) ===')

for variant in VARIANTS:
    blocked = variant == 'within-farm'
    print(f'\\n--- {variant} ---')
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
            y_oof  = cross_val_predict(make_reg(**bp), X, y, cv=cv)
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
            y_eval = np.empty_like(y)
            y_oof  = np.empty_like(y)
            for tr, te in cv.split(X):
                y_tr_c, means = center_y_by_farm(y[tr], farm[tr])
                y_te_c, _     = center_y_by_farm(y[te], farm[te], means)
                model = make_pipeline_reg(**bp).fit(X[tr], y_tr_c)
                y_oof[te]  = model.predict(X[te])
                y_eval[te] = y_te_c
            perm = run_permutation_test(X, center_y_by_farm(y, farm)[0], bp,
                                        key + '__perm', kind='reg', blocked=True)

        rho, _ = spearmanr(y_eval, y_oof)
        REG_RESULTS[(ds, variant)] = dict(
            nested=nested, perm=perm, X_df=X_df, y=y_eval, y_raw=y, y_oof=y_oof,
            spearman=rho, r2=r2_score(y_eval, y_oof),
            mae=mean_absolute_error(y_eval, y_oof), n=len(y), blocked=blocked,
        )
        print(f'  {ds:9s} rho={rho:.3f}  R2={REG_RESULTS[(ds, variant)]["r2"]:.3f}  '
              f'MAE={REG_RESULTS[(ds, variant)]["mae"]:.1f}  p={perm["pvalue"]:.4f}')''')

pathlib.Path(__file__).with_name('cells_part2.json').write_text(
    json.dumps(CELLS, ensure_ascii=False))
print(f'part 2: {len(CELLS)} cells')
