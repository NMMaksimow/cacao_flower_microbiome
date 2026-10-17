# Archived notebook generator

These scripts built `notebooks/03_16S_ITS1_RFC.ipynb` programmatically during the
September 2026 rewrite, when almost every cell changed at once and hand-editing the
notebook JSON would have been error-prone.

**No longer maintained.** The notebook itself is now the single source of truth —
edit `03_16S_ITS1_RFC.ipynb` directly. These files are kept only in case another
wholesale rewrite is ever needed, and they do NOT reflect changes made to the
notebook after 2026-09-13.

`nb_to_script.py` is the exception: it converts the notebook into the batch script
`03_16S_ITS1_RFC.py` used on the HPC, and stays useful.
