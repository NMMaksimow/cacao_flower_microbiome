"""QC for rule group 02: v2 demultiplexed reads must equal the v1 reads.

v1 demultiplexed L1+L2 concatenated, so each v1 file holds the L1 reads followed by
the L2 reads. process_radtags keeps read order and headers, so for every v1 file the
decompressed v1 content must equal the decompressed v2 L1 file followed by the
decompressed v2 L2 file. Contents are compared by md5 (the .gz files themselves may
differ byte-wise). The 16S cmt53 files are also compared with the SRA-deposited
copies, which carry the cmt53/cmtPC correction.

Run by Snakemake (script: directive); paths come from the `snakemake` object.
"""
import glob
import gzip
import hashlib
import os
from concurrent.futures import ProcessPoolExecutor


def md5_of(paths):
    """md5 of the decompressed content of one or more gzip files, read in order."""
    digest = hashlib.md5()
    for path in paths:
        with gzip.open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
    return digest.hexdigest()


def v2_files(name):
    """The L1 and L2 v2 files with this file name (one per lane)."""
    files = []
    for lane in snakemake.params.lanes:
        hits = glob.glob(f"{snakemake.params.demux_dir}/{lane}/*/{name}")
        if len(hits) != 1:
            raise ValueError(f"{name}: expected 1 file in {lane}, found {len(hits)}")
        files.append(hits[0])
    return files


def compare(task):
    source, reference, candidates = task
    return source, reference, md5_of([reference]), md5_of(candidates)


tasks = [
    ("v1", path, v2_files(os.path.basename(path)))
    for path in sorted(glob.glob(f"{snakemake.params.v1_dir}/*/*.fq.gz"))
]
tasks += [
    ("sra", f"{snakemake.params.sra_dir}/{name}", v2_files(name))
    for name in ("cmt53_16S.1.fq.gz", "cmt53_16S.2.fq.gz")
]

with ProcessPoolExecutor(snakemake.threads) as pool:
    results = list(pool.map(compare, tasks))

with open(snakemake.output[0], "w") as out:
    out.write("source\treference_file\tmd5_reference\tmd5_v2_L1_then_L2\tidentical\n")
    for source, reference, md5_ref, md5_v2 in results:
        out.write(f"{source}\t{reference}\t{md5_ref}\t{md5_v2}\t{md5_ref == md5_v2}\n")
