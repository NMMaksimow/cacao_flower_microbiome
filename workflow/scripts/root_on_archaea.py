"""Root the prokaryote 16S tree (track prok) on the archaea, then drop the archaea.

The archaea must form one clade in the unrooted tree. The check: root the tree on
any bacterial tip; the archaea then form one clade exactly when their lowest common
ancestor holds no other tips. The report always lists the number of archaeal and
bacterial ASVs and the result. If the archaea do not form one clade, it also lists
the archaeal ASVs outside the largest archaea-only clade, their taxonomy and the
classes of the nearest bacteria (where they sit), and the script stops without
rooting: a fallback is chosen by hand, not here.

If they form one clade, scikit-bio's root_by_outgroup places the root at the
midpoint of the branch between archaea and bacteria. shear() then keeps only the
bacterial tips, so the root sits at the base of the bacteria.

Usage (rule root_prok_on_archaea, inside the QIIME 2 module for scikit-bio):
    python root_on_archaea.py <unrooted.nwk> <taxonomy.tsv> <report.tsv> \
        <rooted_with_archaea.nwk> <rooted_bacteria.nwk>
"""
import csv
import sys
from collections import Counter

from skbio import TreeNode

tree_path, taxonomy_path, report_path, archaea_out, bacteria_out = sys.argv[1:6]

tree = TreeNode.read(tree_path)
tip_names = {tip.name for tip in tree.tips()}
with open(taxonomy_path) as handle:
    taxonomy = {row["Feature ID"]: row["Taxon"] for row in csv.DictReader(handle, delimiter="\t")
                if row["Feature ID"] in tip_names}
archaea = {name for name, taxon in taxonomy.items() if taxon.startswith("d__Archaea")}
bacteria = tip_names - archaea


def short_taxon(taxon):
    """Phylum and class only, e.g. 'p__Nanobdellota;c__Nanoarchaeia'."""
    return ";".join(part for part in taxon.split(";") if part.startswith(("p__", "c__")))


# Monophyly check on a copy rooted on the branch above one bacterial tip
check = tree.root_at(tree.find(sorted(bacteria)[0]), above=True)
archaea_lca = check.lca(sorted(archaea))
monophyletic = {tip.name for tip in archaea_lca.tips(include_self=True)} == archaea

report = [["archaeal_asvs", len(archaea), "", ""],
          ["bacterial_asvs", len(bacteria), "", ""],
          ["archaea_one_clade", "yes" if monophyletic else "no", "", ""]]

if not monophyletic:
    # Largest clade made of archaea only (counted from the tips up)
    largest = None
    for node in check.postorder():
        if node.is_tip():
            node.n_archaea, node.n_tips = int(node.name in archaea), 1
        else:
            node.n_archaea = sum(child.n_archaea for child in node.children)
            node.n_tips = sum(child.n_tips for child in node.children)
        if node.n_tips == node.n_archaea and (largest is None or node.n_tips > largest.n_tips):
            largest = node
    in_largest = {tip.name for tip in largest.tips(include_self=True)}
    report.append(["largest_archaea_only_clade", len(in_largest), "", ""])
    for name in sorted(archaea - in_largest):
        # Walk up until the clade includes bacteria; report their most common classes
        node = check.find(name).parent
        while not any(tip.name in bacteria for tip in node.tips()):
            node = node.parent
        neighbours = Counter(short_taxon(taxonomy.get(tip.name, "")).split(";")[-1]
                             for tip in node.tips() if tip.name in bacteria)
        report.append(["outside_largest_clade", name, short_taxon(taxonomy[name]),
                       ", ".join(f"{cls} ({n})" for cls, n in neighbours.most_common(3))])

with open(report_path, "w") as out:
    out.write("check\tvalue\ttaxonomy\tnearest_bacteria\n")
    for row in report:
        out.write("\t".join(map(str, row)) + "\n")

if not monophyletic:
    sys.exit(f"The archaea do not form one clade; see {report_path}. No tree was rooted.")

rooted = tree.root_by_outgroup(archaea)
rooted.write(archaea_out)
bacteria_tree = rooted.shear(bacteria)
bacteria_tree.length = None   # the root has no branch above it
bacteria_tree.write(bacteria_out)
