# =============================================================================
# 10. ITS1 classifiers: UNITE v10 (19.02.2025), dynamic, all eukaryotes
# =============================================================================
# Two classifiers, compared on the ITS1 OTUs in QC 07 (as the 16S classifiers in QC 06):
#   s_all          v1's file WITH singletons, trained here exactly as v1 (08a):
#                  import FASTA and taxonomy, fit-classifier-naive-bayes with
#                  defaults, full-length references (no region extraction)
#   no_singletons  community classifier for QIIME 2 2026.4 (unite-train), downloaded
# Both go to the shared databases folder (classifiers/qiime2-2026.4_sklearn-1.7.1/).

UNITE_REF_DIR = f"{DB_DIR}/{config['unite']['reference_dir']}"
UNITE_CLF_DIR = f"{DB_DIR}/{config['unite']['classifier_dir']}"
UNITE = {name: f"{UNITE_CLF_DIR}/{file}" for name, file in config["unite"]["classifiers"].items()}

localrules: download_unite_no_singletons


rule train_unite_s_all:
    input:
        fasta=f"{UNITE_REF_DIR}/{config['unite']['fasta']}",
        taxonomy=f"{UNITE_REF_DIR}/{config['unite']['taxonomy']}",
    output:
        UNITE["s_all"],
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
    threads: 1
    resources:
        partition="bigmem",   # 1 CPU with 240 GB fits easily on a 2 TB node
        mem_mb=240000,
        runtime=1440,
        tmp_gb=80,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        work=$(mktemp -d)
        {params.r_isolation} qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {input.fasta} \
            --output-path $work/sequences.qza
        {params.r_isolation} qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {input.taxonomy} \
            --output-path $work/taxonomy.qza
        {params.r_isolation} qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads $work/sequences.qza \
            --i-reference-taxonomy $work/taxonomy.qza \
            --o-classifier {output} \
            --verbose
        rm -r $work
        """


rule download_unite_no_singletons:
    output:
        UNITE["no_singletons"],
    params:
        url=config["unite"]["community_url"],
    shell:
        r"""
        set -euo pipefail
        curl -sSfL -o {output} {params.url}
        md5sum {output} | sed 's|  .*/|  |' > {output}.md5
        """

