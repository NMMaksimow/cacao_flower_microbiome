# =============================================================================
# 10. ITS1 classifiers: UNITE v10 (19.02.2025), dynamic, all eukaryotes
# =============================================================================
# Two classifiers, compared in 2.7 (F-score on a common test set + on-data checks):
#   s_all          v1's file WITH singletons, trained here exactly as v1 (08a):
#                  import FASTA and taxonomy, fit-classifier-naive-bayes with
#                  defaults, full-length references (no region extraction)
#   no_singletons  community classifier for QIIME 2 2026.4 (unite-train), downloaded
# Both go to the shared databases folder (classifiers/qiime2-2026.4_sklearn-1.7.1/).

UNITE_REF_DIR = f"{DB_DIR}/{config['unite']['reference_dir']}"
UNITE_CLF_DIR = f"{DB_DIR}/{config['unite']['classifier_dir']}"
UNITE = {name: f"{UNITE_CLF_DIR}/{file}" for name, file in config["unite"]["classifiers"].items()}

localrules: download_unite_no_singletons, unite_test_set


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
        mem_mb=240000,
        runtime=2880,
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


rule unite_test_set:
    input:
        fasta=f"{UNITE_REF_DIR}/{config['unite']['fasta']}",
        taxonomy=f"{UNITE_REF_DIR}/{config['unite']['taxonomy']}",
    output:
        fasta=f"{OUT_DIR}/taxonomy/unite_test_set/test_sequences.fasta",
        taxonomy=f"{OUT_DIR}/taxonomy/unite_test_set/test_taxonomy.tsv",
    params:
        script=f"{workflow.basedir}/scripts/unite_test_set.py",
        n=config["unite"]["test_set_size"],
        seed=config["dada2"]["seed"],
    shell:
        r"""
        python3 {params.script} {input.fasta} {input.taxonomy} {params.n} {params.seed} \
            {output.fasta} {output.taxonomy}
        """
