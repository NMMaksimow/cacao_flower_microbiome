"""Convert 03_16S_ITS1_RFC.ipynb into a batch-friendly .py script.

Done by hand rather than with `jupyter nbconvert` so that jupyter is not needed on
the cluster and so the two batch-specific changes are explicit:
  * a headless matplotlib backend, chosen before pyplot is imported;
  * QUICK driven by an environment variable, so a cheap smoke run can verify the
    environment before a multi-hour job is submitted.
"""
import json, pathlib, re

NB  = pathlib.Path('/mnt/c/Users/nmmak/Desktop/cacao_flower_microbiome/notebooks/03_16S_ITS1_RFC.ipynb')
OUT = pathlib.Path('/mnt/c/Users/nmmak/Desktop/cacao_flower_microbiome/notebooks/03_16S_ITS1_RFC.py')

HEADER = '''#!/usr/bin/env python3
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

'''

cells = [c for c in json.loads(NB.read_text())['cells'] if c['cell_type'] == 'code']
parts = [HEADER]

for i, cell in enumerate(cells):
    src = ''.join(cell['source'])
    if i == 0:
        # Let the batch job pick the profile without editing the file.
        src = re.sub(r'^QUICK\s*=\s*(True|False).*$',
                     "QUICK = os.environ.get('RFC_QUICK', '0') == '1'",
                     src, count=1, flags=re.M)
        assert "os.environ.get('RFC_QUICK'" in src, 'QUICK line not found'
    parts.append(f'# ===== cell {i} ' + '=' * 58 + f'\n{src}\n')

OUT.write_text('\n'.join(parts))
print(f'wrote {OUT}  ({OUT.stat().st_size / 1024:.0f} KB, {len(cells)} cells)')

import ast
ast.parse(OUT.read_text())
print('syntax OK')
