#!/bin/bash

echo "=== DEMULTIPLEXING COMPARISON ==="
echo "Comparing original (-c -q -r) vs new (no quality filters)"
echo ""

# Function to count reads in a directory
count_reads() {
    local dir=$1
    local label=$2
    echo "Processing $label..."
    
    total_reads=0
    file_count=0
    
    for file in "$dir"/*.fq.gz; do
        if [[ -f "$file" ]]; then
            reads=$(zcat "$file" | wc -l)
            reads=$((reads / 4))  # Convert lines to read count
            total_reads=$((total_reads + reads))
            file_count=$((file_count + 1))
        fi
    done
    
    echo "$label - Files: $file_count, Total reads: $total_reads"
    if [[ $file_count -gt 0 ]]; then
        avg_reads=$((total_reads / file_count))
        echo "$label - Average reads per file: $avg_reads"
    fi
    echo ""
}

# Compare 16S
echo "=== 16S COMPARISON ==="
count_reads "qiime2/import/demux/demultiplexed_sample_files/16S" "Original 16S"
count_reads "qiime2/import/demux/demultiplexed_sample_files_no_cqr/16S" "New 16S"

# Compare ITS1
echo "=== ITS1 COMPARISON ==="
count_reads "qiime2/import/demux/demultiplexed_sample_files/ITS1" "Original ITS1"
count_reads "qiime2/import/demux/demultiplexed_sample_files_no_cqr/ITS1" "New ITS1"
