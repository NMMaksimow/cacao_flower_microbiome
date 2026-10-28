# =============================================================================
# 13. 16S phylogeny: three trees ("tracks") that differ only in input set and root
# =============================================================================
# Tree method as v1 (qiime2/scripts/12): align-to-tree-mafft-fasttree with defaults
# (MAFFT, masking, FastTree, midpoint root). Input sets (config phylogeny_tracks):
#   all   every ASV (v1 method on v2 data), midpoint root
#   prok  bacteria + archaea without organelles, rooted on the archaea (outgroup)
#   bact  bacteria without organelles (= scripts/01a filter), midpoint root
# bact is the primary tree; all and prok are sensitivity analyses (decision rules
# in docs/v2_rerun_plan.md). ITS1 has no phylogeny.
#
# q2-phylogeny can only midpoint-root, so the prok tree is rooted by
# scripts/root_on_archaea.py (scikit-bio root_by_outgroup). It stops, with a
# report, if the archaea do not form one clade.

PHYLO_TRACKS = ["all"] + list(config["phylogeny_tracks"])
MAIN_TAXONOMY_16S = f"{OUT_DIR}/taxonomy/CFM_16S_taxonomy_{config['taxonomy_main']['16S']}"


rule filter_16S_for_tree:
    input:
        reads=f"{OUT_DIR}/denoise/CFM_16S_dada2_repseqs.qza",
        taxonomy=f"{MAIN_TAXONOMY_16S}.qza",
    output:
        f"{OUT_DIR}/phylogeny/{{track}}/CFM_16S_repseqs.qza",
    wildcard_constraints:
        track="prok|bact",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
        include=lambda wildcards: config["phylogeny_tracks"][wildcards.track]["include"],
        exclude=lambda wildcards: config["phylogeny_tracks"][wildcards.track]["exclude"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime taxa filter-seqs \
            --i-sequences {input.reads} \
            --i-taxonomy {input.taxonomy} \
            --p-include '{params.include}' \
            --p-exclude '{params.exclude}' \
            --p-mode contains \
            --o-filtered-sequences {output}
        """


rule phylogeny_16S:
    input:
        lambda wildcards: (f"{OUT_DIR}/denoise/CFM_16S_dada2_repseqs.qza" if wildcards.track == "all"
                           else f"{OUT_DIR}/phylogeny/{wildcards.track}/CFM_16S_repseqs.qza"),
    output:
        alignment=f"{OUT_DIR}/phylogeny/{{track}}/CFM_16S_aligned_repseqs.qza",
        masked=f"{OUT_DIR}/phylogeny/{{track}}/CFM_16S_masked_aligned_repseqs.qza",
        unrooted=f"{OUT_DIR}/phylogeny/{{track}}/CFM_16S_unrooted_tree.qza",
        rooted=f"{OUT_DIR}/phylogeny/{{track}}/CFM_16S_rooted_tree.qza",
    wildcard_constraints:
        track="all|prok|bact",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: config["phylogeny_threads"]
    resources:
        mem_mb=32000,
        runtime=360,   # v1: 58 min for ~29k ASVs on 16 cores
        tmp_gb=20,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input} \
            --p-n-threads {threads} \
            --o-alignment {output.alignment} \
            --o-masked-alignment {output.masked} \
            --o-tree {output.unrooted} \
            --o-rooted-tree {output.rooted}
        """


rule root_prok_on_archaea:
    input:
        tree=f"{OUT_DIR}/phylogeny/prok/CFM_16S_unrooted_tree.qza",
        taxonomy=f"{MAIN_TAXONOMY_16S}.tsv",
    output:
        archaea=f"{OUT_DIR}/phylogeny/prok/CFM_16S_rooted_tree_archaea.nwk",
        bacteria=f"{OUT_DIR}/phylogeny/prok/CFM_16S_rooted_tree_bacteria.nwk",
    # The report is a log: Snakemake keeps logs of failed jobs, so it survives
    # when the script stops because the archaea are not one clade.
    log:
        f"{OUT_DIR}/qc/08_tree_archaea_outgroup.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
        script=f"{workflow.basedir}/scripts/root_on_archaea.py",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        export_dir=$(mktemp -d)
        {params.r_isolation} qiime tools export --input-path {input.tree} --output-path $export_dir
        python {params.script} $export_dir/tree.nwk {input.taxonomy} {log} \
            {output.archaea} {output.bacteria}
        rm -r $export_dir
        """
