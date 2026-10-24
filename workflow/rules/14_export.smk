# =============================================================================
# 14. Exports for R: plain files in qiime2_v2/export/, v1 folder and file names
# =============================================================================
# Folder names as v1's qiime2/export/, so the R scripts switch mainly by path:
#   CFM_<marker>_dada2_table/feature-table.biom, ..._repseqs/dna-sequences.fasta
#   CFM_ITS1_otu97_{table,repseqs}/...
#   CFM_16S_rooted_tree_{bact,all,prok}/tree.nwk   (bact = primary; v1 had one tree)
#   CFM_16S_taxonomy_<classifier>/taxonomy.tsv     (v1: ..._diverse_weighted)
#   CFM_ITS1_taxonomy/, CFM_ITS1_otu97_taxonomy/   (main ITS1 classifier)
#   RUN_INFO.txt                                    (provenance of all of the above)
# The taxonomy of the main classifiers (config taxonomy_main) is copied from
# qiime2_v2/taxonomy/. The ITS1 taxonomy and RUN_INFO wait until ITS1 is chosen.

EXPORT_DIR = f"{OUT_DIR}/export"
EXPORT_FEATURES = [
    f"{EXPORT_DIR}/CFM_{name}_table/feature-table.biom"
    for name in ["16S_dada2", "ITS1_dada2", "ITS1_otu97"]
] + [
    f"{EXPORT_DIR}/CFM_{name}_repseqs/dna-sequences.fasta"
    for name in ["16S_dada2", "ITS1_dada2", "ITS1_otu97"]
]
EXPORT_TREES = expand(f"{EXPORT_DIR}/CFM_16S_rooted_tree_{{track}}/tree.nwk", track=PHYLO_TRACKS)

MAIN_16S = config["taxonomy_main"]["16S"]
MAIN_ITS1 = config["taxonomy_main"]["ITS1"]
TAXONOMY_EXPORT = {f"CFM_16S_taxonomy_{MAIN_16S}": f"{MAIN_TAXONOMY_16S}.tsv"}
if MAIN_ITS1:
    TAXONOMY_EXPORT["CFM_ITS1_taxonomy"] = f"{OUT_DIR}/taxonomy/CFM_ITS1_asv_taxonomy_{MAIN_ITS1}.tsv"
    TAXONOMY_EXPORT["CFM_ITS1_otu97_taxonomy"] = f"{OUT_DIR}/taxonomy/CFM_ITS1_otu97_taxonomy_{MAIN_ITS1}.tsv"
EXPORT_TAXONOMY = [f"{EXPORT_DIR}/{name}/taxonomy.tsv" for name in TAXONOMY_EXPORT]

localrules: export_tree_prok, export_taxonomy, run_info


rule export_artifact:
    input:
        f"{OUT_DIR}/denoise/{{name}}.qza",
    output:
        f"{EXPORT_DIR}/{{name}}/{{file}}",
    wildcard_constraints:
        name="CFM_[^/]+",
        file="feature-table.biom|dna-sequences.fasta",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime tools export --input-path {input} --output-path $(dirname {output})
        """


rule export_tree:
    input:
        f"{OUT_DIR}/phylogeny/{{track}}/CFM_16S_rooted_tree.qza",
    output:
        f"{EXPORT_DIR}/CFM_16S_rooted_tree_{{track}}/tree.nwk",
    wildcard_constraints:
        track="all|bact",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime tools export --input-path {input} --output-path $(dirname {output})
        """


# prok: the archaea-rooted tree with the archaea dropped (rule root_prok_on_archaea)
rule export_tree_prok:
    input:
        f"{OUT_DIR}/phylogeny/prok/CFM_16S_rooted_tree_bacteria.nwk",
    output:
        f"{EXPORT_DIR}/CFM_16S_rooted_tree_prok/tree.nwk",
    shell:
        "cp {input} {output}"


rule export_taxonomy:
    input:
        lambda wildcards: TAXONOMY_EXPORT[wildcards.name],
    output:
        f"{EXPORT_DIR}/{{name}}/taxonomy.tsv",
    wildcard_constraints:
        name="CFM_[^/]+",
    shell:
        "cp {input} {output}"


rule run_info:
    input:
        exports=EXPORT_FEATURES + EXPORT_TREES + EXPORT_TAXONOMY,
        versions=f"{OUT_DIR}/00_env/tool_versions.txt",
    output:
        f"{EXPORT_DIR}/RUN_INFO.txt",
    params:
        code_dir=workflow.basedir,
        qiime2=config["module_qiime2"],
        main_16S=MAIN_16S,
        main_ITS1=MAIN_ITS1,
        silva_main=lambda wildcards: SILVA144[MAIN_16S],
        unite_main=lambda wildcards: UNITE[MAIN_ITS1],
        unite_source=lambda wildcards: (config["unite"]["community_url"] if MAIN_ITS1 == "no_singletons"
                                        else "trained in this pipeline (rule train_unite_s_all)"),
        unite_doi=config["unite"]["doi"],
        tracks=" ".join(PHYLO_TRACKS),
    shell:
        r"""
        {{
        echo "CFM v2 exports (qiime2_v2/export/)"
        echo "date: $(date -Iseconds)"
        echo "code: $(git -C {params.code_dir} rev-parse HEAD) (branch $(git -C {params.code_dir} rev-parse --abbrev-ref HEAD))"
        echo "config: {params.code_dir}/config_CFM_v2.yaml"
        echo "QIIME 2 module: {params.qiime2}"
        echo
        echo "16S classifier: SILVA 144 V4 515F/806R, {params.main_16S}"
        echo "  file: {params.silva_main}"
        echo "  DOI: $(cat {params.silva_main}.doi)"
        echo "  md5: $(md5sum {params.silva_main} | cut -d' ' -f1)"
        echo "ITS1 classifier: UNITE v10 19.02.2025, dynamic, all eukaryotes, {params.main_ITS1}"
        echo "  file: {params.unite_main}"
        echo "  source: {params.unite_source}"
        echo "  UNITE DOI: {params.unite_doi}"
        echo "  md5: $(md5sum {params.unite_main} | cut -d' ' -f1)"
        echo "16S trees: {params.tracks} (primary: bact)"
        echo
        echo "Files (md5):"
        (cd $(dirname {output}) && md5sum */*)
        echo
        echo "Tool versions (qiime2_v2/00_env/tool_versions.txt):"
        cat {input.versions}
        }} > {output}
        """
