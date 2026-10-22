"""Common test set for comparing the two UNITE classifiers (F-score per rank).

Draws a fixed-seed random sample of non-singleton reference sequences (header suffix
not "_singleton") from the UNITE s_all FASTA. Both classifiers (s_all and
no_singletons) were trained on these sequences, so the test measures how much the
singletons disturb the assignment of well-known taxa, not accuracy on novel taxa.
The truth is the UNITE taxonomy with the 8th level (sh__) removed, because the
no-singleton classifier was trained without it.

Usage (rule unite_test_set):
    python unite_test_set.py <unite.fasta> <unite_taxonomy.txt> <n> <seed> \
        <test.fasta> <test_taxonomy.tsv>
"""
import random
import sys

fasta_in, taxonomy_in, n, seed, fasta_out, taxonomy_out = sys.argv[1:7]

sequences, name = {}, None
for line in open(fasta_in):
    if line.startswith(">"):
        name = line[1:].strip()
        sequences[name] = []
    else:
        sequences[name].append(line.strip())

candidates = sorted(i for i in sequences if not i.endswith("_singleton"))
chosen = sorted(random.Random(int(seed)).sample(candidates, int(n)))

taxonomy = {}
for line in open(taxonomy_in):
    feature_id, taxon = line.rstrip("\n").split("\t")[:2]
    taxonomy[feature_id] = ";".join(r for r in taxon.split(";") if not r.startswith("sh__"))

with open(fasta_out, "w") as out:
    for i in chosen:
        out.write(f">{i}\n{''.join(sequences[i])}\n")
with open(taxonomy_out, "w") as out:
    out.write("Feature ID\tTaxon\n")
    for i in chosen:
        out.write(f"{i}\t{taxonomy[i]}\n")
print(f"{len(candidates)} non-singleton references, {len(chosen)} drawn")
