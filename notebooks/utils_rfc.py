"""Shared input for the random forest scripts (notebooks 03, 03b, 03c).

The CLR features are not computed in Python. They come from
scripts/08c_16S_ITS1_clr_export.R, which applies the same function as 08a/08b
(clr_czm() in scripts/utils_clr.R: prevalence filter, zCompositions::cmultRepl
with the CZM method, CLR). One computation, so the PCoA, the PERMANOVA and the
random forest all see the same CLR values.
"""

from pathlib import Path

import pandas as pd

# Sample set -> tag in the 08c file names (SET_TAG in scripts/utils_paths.R)
SET_TAG = {'main': '_288', 'with_deviations': '_294'}

# 08c prevalence configurations, see the header of 08c
CLR_CONFIGS = ['all10', 'all05', 'all20', 'unbagged05', 'unbagged10']


def load_clr(marker, results_dir, sample_set='main', config='all10', counts=None):
    """Read one 08c CLR matrix. Returns a DataFrame, samples x features.

    marker      '16S' or 'ITS1'
    results_dir results/<version>
    config      prevalence configuration of 08c; 'all10' is the main one
    counts      optional count table from the biom (features x samples); when given,
                the CLR features must be a subset of its features and the samples
                must be the same, otherwise the R and Python inputs have drifted apart
    """
    path = (Path(results_dir) / 'tables' / 'clr' /
            f'08c_{marker}_clr_czm_{config}{SET_TAG[sample_set]}.csv.gz')
    if not path.exists():
        raise FileNotFoundError(f'{path} -- run scripts/08c_16S_ITS1_clr_export.R '
                                f'for sample set "{sample_set}" first')
    X = pd.read_csv(path, index_col='sample_id', dtype={'sample_id': str})
    X.columns = X.columns.astype(str)

    if counts is not None:
        extra = X.columns.difference(counts.index)
        if len(extra):
            raise ValueError(f'{path.name}: {len(extra)} features not in the biom table')
        if set(X.index) != set(counts.columns):
            raise ValueError(f'{path.name}: samples differ from the biom table '
                             f'({len(X)} vs {counts.shape[1]})')
    return X
