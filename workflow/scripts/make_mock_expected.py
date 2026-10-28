"""Expected amplicon sequences of the ZymoBIOMICS D6305 mock community.

In-silico PCR on the Zymo reference archive (ZymoBIOMICS.STD.refseq.v3.zip):
  16S V4 (515F/806R) on every 16S copy in ssrRNAs/ (8 bacteria);
  ITS1 (ITS1F/ITS2) on the complete genomes of the two yeasts.
A primer site may differ from the primer by up to 2 mismatches, as cutadapt allows
for these primers (error rate 0.1). The output sequences are the inserts between
the primers, i.e. what v2 reads look like after primer trimming. Identical copies
are collapsed; the header gives the number of copies.

Run once on the laptop; the FASTA files are committed in workflow/resources/mock/:
    python3 workflow/scripts/make_mock_expected.py \
        data/ZymoBIOMICS.STD.refseq.v3.zip workflow/resources/mock
"""
import io
import re
import sys
import zipfile
from collections import Counter

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
         "V": "ACG", "N": "ACGT"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

# Primer pairs: forward primer and reverse complement of the reverse primer
PRIMERS = {
    "16S_V4": ("GTGYCAGCMGCCGCGGTAA", "ATTAGAWACCCBNGTAGTCC"),   # 515F / RC(806R)
    "ITS1": ("CTTGGTCATTTAGAGGAAGTAA", "GCATCGATGAAGAACGCAGC"),  # ITS1F / RC(ITS2)
}
MAX_MISMATCHES = 2
MAX_INSERT = 1000


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(handle):
    """Yield (name, sequence) from a FASTA file opened in binary mode."""
    name, chunks = None, []
    for line in io.TextIOWrapper(handle):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip().upper())
    if name is not None:
        yield name, "".join(chunks)


def mismatches(primer, site):
    return sum(base not in IUPAC[p] for p, base in zip(primer, site))


def primer_sites(seq, primer):
    """Start positions where the primer matches with <= MAX_MISMATCHES.

    With at most 2 mismatches, at least one third of the primer matches exactly,
    so exact matches of the thirds are used as seeds before the full comparison.
    """
    third = len(primer) // 3
    seeds = [(0, primer[:third]), (third, primer[third:2 * third]), (2 * third, primer[2 * third:])]
    starts = set()
    for offset, part in seeds:
        pattern = re.compile("".join(f"[{IUPAC[b]}]" for b in part))
        for match in pattern.finditer(seq):
            start = match.start() - offset
            if 0 <= start <= len(seq) - len(primer):
                if mismatches(primer, seq[start:start + len(primer)]) <= MAX_MISMATCHES:
                    starts.add(start)
    return sorted(starts)


def inserts(seq, fwd, rev_rc):
    """Inserts between each forward site and the nearest downstream reverse site."""
    found = []
    for strand in (seq, reverse_complement(seq)):
        rev_sites = primer_sites(strand, rev_rc)
        for f in primer_sites(strand, fwd):
            start = f + len(fwd)
            downstream = [r for r in rev_sites if start < r <= start + MAX_INSERT]
            if downstream:
                found.append(strand[start:downstream[0]])
    return found


def organism(path):
    """Genus_species from a file name such as ssrRNAs/Bacillus_subtilis_16S_170923.fasta."""
    return "_".join(path.split("/")[-1].split("_")[:2])


def write_fasta(path, records):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n{seq}\n")


archive_path, out_dir = sys.argv[1], sys.argv[2]
archive = zipfile.ZipFile(archive_path)
sources = {
    "16S_V4": [n for n in archive.namelist() if "/ssrRNAs/" in n and "_16S_" in n],
    "ITS1": [n for n in archive.namelist() if "/Genomes/" in n
             and re.search("Saccharomyces|Cryptococcus", n)],
}
for amplicon, files in sources.items():
    fwd, rev_rc = PRIMERS[amplicon]
    records = []
    for path in sorted(files):
        copies = Counter(insert for _, seq in read_fasta(archive.open(path))
                         for insert in inserts(seq, fwd, rev_rc))
        for number, (seq, count) in enumerate(copies.most_common(), start=1):
            records.append((f"{organism(path)}|variant{number}|copies={count}|length={len(seq)}", seq))
        print(f"{amplicon}\t{organism(path)}\tcopies={sum(copies.values())}\tvariants={len(copies)}")
    write_fasta(f"{out_dir}/D6305_v3_{amplicon}.fasta", records)
