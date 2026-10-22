# =============================================================================
# 09. 16S taxonomy with the four official SILVA 144 classifiers
# =============================================================================
# The classifiers are downloaded once into the shared databases folder, each with
# its .md5 (checked) and .doi file. Downloads run on the login node ("localrules"),
# because compute nodes may have no internet access. classify-sklearn uses the
# default confidence threshold 0.7, as v1.

SILVA144_DIR = f"{DB_DIR}/{config['silva144']['dir']}"
SILVA144 = {name: f"{SILVA144_DIR}/{os.path.basename(path)}"
            for name, path in config["silva144"]["classifiers"].items()}

localrules: download_silva_classifiers


rule download_silva_classifiers:
    output:
        list(SILVA144.values()),
    params:
        base_url=config["silva144"]["base_url"],
        paths=" ".join(config["silva144"]["classifiers"].values()),
        out_dir=SILVA144_DIR,
    shell:
        r"""
        set -euo pipefail
        cd {params.out_dir}
        for path in {params.paths}; do
            for suffix in "" .md5 .doi; do
                curl -sSfL -o "$(basename $path)$suffix" "{params.base_url}/$path$suffix"
            done
            md5sum -c "$(basename $path).md5"
        done
        """


rule classify_16S:
    input:
        classifier=lambda wildcards: SILVA144[wildcards.classifier],
        reads=f"{OUT_DIR}/denoise/CFM_16S_dada2_repseqs.qza",
    output:
        qza=f"{OUT_DIR}/taxonomy/CFM_16S_taxonomy_{{classifier}}.qza",
        tsv=f"{OUT_DIR}/taxonomy/CFM_16S_taxonomy_{{classifier}}.tsv",
    params:
        base=config["module_base"],
        qiime2=config["module_qiime2"],
        r_isolation=config["r_isolation"],
        tmp_dir=config["tmp_dir"],
        n_jobs=config["classify_n_jobs"],
    threads: config["classify_n_jobs"]
    resources:
        mem_mb=96000,
        runtime=360,
        tmp_gb=20,
    shell:
        r"""
        set -euo pipefail
        module load {params.base} {params.qiime2}
        export TMPDIR={params.tmp_dir}
        {params.r_isolation} qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.reads} \
            --p-n-jobs {params.n_jobs} \
            --o-classification {output.qza}
        export_dir=$(mktemp -d)
        {params.r_isolation} qiime tools export --input-path {output.qza} --output-path $export_dir
        mv $export_dir/taxonomy.tsv {output.tsv}
        rm -r $export_dir
        """
