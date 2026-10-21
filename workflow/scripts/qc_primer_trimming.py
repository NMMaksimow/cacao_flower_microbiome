"""QC for rule group 03: summary of primer trimming per lane and marker.

Reads the cutadapt JSON report of every sample and writes two tables:
  1. per lane x marker (v2) and per marker (v1, from logs/05_primer_trimming.out):
     read pairs in and written, pairs discarded for a missing primer or for being
     too short, and how many reads had the 5' primer and the 3' read-through
     (reverse-complement primer) removed;
  2. per sample and marker: pairs in and written per lane, and the retention.

Run by Snakemake (script: directive); paths come from the `snakemake` object.
"""
import json
import re
from collections import defaultdict

COLUMNS = ["pairs_in", "pairs_written", "no_primer", "too_short",
           "r1_5prime", "r1_3prime", "r2_5prime", "r2_3prime"]


def v2_counts(path):
    """Counts of one sample from its cutadapt JSON report (linked adapters)."""
    report = json.load(open(path))
    reads, r1, r2 = report["read_counts"], report["adapters_read1"][0], report["adapters_read2"][0]
    return {
        "pairs_in": reads["input"],
        "pairs_written": reads["output"],
        "no_primer": reads["filtered"]["discard_untrimmed"],
        "too_short": reads["filtered"]["too_short"],
        "r1_5prime": r1["five_prime_end"]["matches"],
        "r1_3prime": r1["three_prime_end"]["matches"],
        "r2_5prime": r2["five_prime_end"]["matches"],
        "r2_3prime": r2["three_prime_end"]["matches"],
    }


def v1_counts(log_path, primers):
    """Per-marker totals from the v1 log: pairs in/written and 5' primer trims."""
    totals = {marker: defaultdict(int) for marker in primers}
    marker = None
    for line in open(log_path):
        if line.startswith("Command line parameters:"):
            marker = re.search(r"demultiplexed_sample_files/(\w+)/", line).group(1)
        elif line.startswith("Total read pairs processed:"):
            totals[marker]["pairs_in"] += int(line.split()[-1].replace(",", ""))
        elif line.startswith("Pairs written (passing filters):"):
            totals[marker]["pairs_written"] += int(line.split()[4].replace(",", ""))
        elif "Type: regular 5'" in line:
            sequence = line.split(";")[0].split()[-1]
            trimmed = int(re.search(r"Trimmed: (\d+) times", line).group(1))
            read = "r1_5prime" if sequence == primers[marker]["fwd"] else "r2_5prime"
            totals[marker][read] += trimmed
    return totals


per_group = defaultdict(lambda: defaultdict(int))
per_sample = defaultdict(dict)
for path in snakemake.input:
    lane, marker, filename = path.split("/")[-3:]
    sample = filename.replace(".cutadapt.json", "")
    counts = v2_counts(path)
    for column in COLUMNS:
        per_group[("v2", lane, marker)][column] += counts[column]
    per_sample[(sample, marker)][lane] = counts

v1 = v1_counts(snakemake.params.v1_log, snakemake.params.primers)
for marker, counts in v1.items():
    per_group[("v1", "L1+L2", marker)] = counts

with open(snakemake.output.summary, "w") as out:
    out.write("run\tlane\tmarker\t" + "\t".join(COLUMNS) + "\tpct_written\tpct_r1_readthrough\n")
    for (run, lane, marker), counts in sorted(per_group.items(), key=lambda x: (x[0][2], x[0][0], x[0][1])):
        values = [counts.get(column, "") for column in COLUMNS]
        pct_written = 100 * counts["pairs_written"] / counts["pairs_in"]
        readthrough = 100 * counts["r1_3prime"] / counts["pairs_in"] if "r1_3prime" in counts else ""
        readthrough = f"{readthrough:.2f}" if readthrough != "" else ""
        out.write(f"{run}\t{lane}\t{marker}\t" + "\t".join(map(str, values))
                  + f"\t{pct_written:.2f}\t{readthrough}\n")

with open(snakemake.output.samples, "w") as out:
    out.write("sample\tmarker\tL1_in\tL1_written\tL2_in\tL2_written\tpct_written\n")
    for (sample, marker), lanes in sorted(per_sample.items()):
        pairs_in = sum(lanes[lane]["pairs_in"] for lane in lanes)
        written = sum(lanes[lane]["pairs_written"] for lane in lanes)
        pct = f"{100 * written / pairs_in:.2f}" if pairs_in else "NA"
        out.write(f"{sample}\t{marker}\t{lanes['L1']['pairs_in']}\t{lanes['L1']['pairs_written']}\t"
                  f"{lanes['L2']['pairs_in']}\t{lanes['L2']['pairs_written']}\t{pct}\n")
