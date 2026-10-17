"""Assemble the four cell parts into notebooks/03_16S_ITS1_RFC.ipynb."""
import json, pathlib, sys

HERE = pathlib.Path(__file__).parent
OUT  = pathlib.Path('/mnt/c/Users/nmmak/Desktop/cacao_flower_microbiome/notebooks/03_16S_ITS1_RFC.ipynb')

cells_spec = []
for part in (1, 2, 3, 4):
    p = HERE / f'cells_part{part}.json'
    if not p.exists():
        sys.exit(f'missing {p} - run build_nb_part{part}.py first')
    cells_spec.extend(json.loads(p.read_text()))

# Preserve the kernel/language metadata of the existing notebook so the file keeps
# pointing at the same conda environment the user already has set up.
old = json.loads(OUT.read_text()) if OUT.exists() else {}
metadata = old.get('metadata', {})

def to_source(src):
    """nbformat stores source as a list of lines, each keeping its trailing newline."""
    lines = src.split('\n')
    return [l + '\n' for l in lines[:-1]] + [lines[-1]]

cells = []
for kind, src in cells_spec:
    cell = {'cell_type': kind, 'metadata': {}, 'source': to_source(src)}
    if kind == 'code':
        cell['execution_count'] = None
        cell['outputs'] = []
    cells.append(cell)

nb = {'cells': cells, 'metadata': metadata,
      'nbformat': old.get('nbformat', 4), 'nbformat_minor': old.get('nbformat_minor', 5)}

OUT.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + '\n')
print(f'wrote {OUT}')
print(f'cells: {len(cells)}  ({sum(c["cell_type"]=="code" for c in cells)} code, '
      f'{sum(c["cell_type"]=="markdown" for c in cells)} markdown)')
print(f'kernel: {metadata.get("kernelspec", {}).get("name", "(none carried over)")}')
