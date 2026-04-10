#!/bin/bash
#SBATCH --job-name=download_wemheuer
#SBATCH -p shared
#SBATCH --mem=2G
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=logs/01_download_seqs_wemheuer.out
#SBATCH --error=logs/01_download_seqs_wemheuer.err

# Download PRJNA594470 (Wemheuer et al. 2020 — cacao leaf endophytes, Cameroon)
# 128 runs × 2 files (R1 + R2) = 256 FASTQ.gz, ~1.6 GB total
#
# Run from project root:
#   sbatch qiime2_wemheuer/scripts/01_download_seqs_wemheuer.sh
#
# Input:  data_wemheuer/metadata.tsv
# Output: data_wemheuer/

echo "=== Download Wemheuer PRJNA594470 ==="
echo "Date:     $(date)"
echo "Host:     $(hostname)"
echo "Job ID:   ${SLURM_JOB_ID}"
echo ""

# ── Paths (relative to project root = SLURM_SUBMIT_DIR) ──────────────────────
METADATA="data_wemheuer/metadata.tsv"
OUTDIR="data_wemheuer"
URL_LIST="logs/01_download_seqs_wemheuer_urls.txt"

mkdir -p "${OUTDIR}" logs

# ── Sanity check ──────────────────────────────────────────────────────────────
if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: ${METADATA} not found."
    echo "Expected at: data_wemheuer/metadata.tsv"
    exit 1
fi

# ── Extract ENA HTTP URLs by column header ────────────────────────────────────
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "ena_fastq_http_1") c1 = i
        if ($i == "ena_fastq_http_2") c2 = i
    }
    next
}
{
    if (c1 && $c1 != "" && $c1 != "-") print $c1
    if (c2 && $c2 != "" && $c2 != "-") print $c2
}
' "${METADATA}" > "${URL_LIST}"

N_URLS=$(wc -l < "${URL_LIST}")
echo "URLs to download: ${N_URLS}"

if [[ "${N_URLS}" -eq 0 ]]; then
    echo "ERROR: No URLs found. Check ${METADATA} column names."
    exit 1
fi

# ── Parallel wget (8 jobs) ────────────────────────────────────────────────────
module load parallel 2>/dev/null || true

echo "Downloading to: ${OUTDIR}/"
echo ""

parallel --jobs 8 --bar \
    wget -q -c --tries=5 --waitretry=2 -P "${OUTDIR}" {} \
    :::: "${URL_LIST}"

# ── Verify ────────────────────────────────────────────────────────────────────
N_DONE=$(find "${OUTDIR}" -name "*.fastq.gz" | wc -l)
echo ""
echo "=== Summary ==="
echo "Expected : ${N_URLS}"
echo "Found    : ${N_DONE}"
echo "Date     : $(date)"

if [[ "${N_DONE}" -lt "${N_URLS}" ]]; then
    echo "WARNING: $((N_URLS - N_DONE)) file(s) missing — check logs/01_download_seqs_wemheuer.err"
    exit 1
fi

echo "All files downloaded successfully."
