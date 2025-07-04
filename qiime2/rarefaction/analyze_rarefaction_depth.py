#!/usr/bin/env python3
"""
Script to extract key information from QIIME2 feature table summaries
to help determine optimal rarefaction depths.
"""

import pandas as pd
import zipfile
import json
import sys
import os

def extract_sample_stats(qzv_file):
    """Extract sample frequency statistics from QIIME2 visualization file."""
    with zipfile.ZipFile(qzv_file, 'r') as zip_ref:
        # Find the data file
        file_list = zip_ref.namelist()
        data_file = None
        for file in file_list:
            if file.endswith('sample-frequency-detail.csv'):
                data_file = file
                break
        
        if data_file:
            with zip_ref.open(data_file) as f:
                df = pd.read_csv(f)
                return df
        else:
            print(f"Could not find sample frequency data in {qzv_file}")
            return None

def analyze_rarefaction_depth(df, dataset_name):
    """Analyze sample depths and suggest rarefaction parameters."""
    if df is None:
        return
    
    print(f"\n=== {dataset_name} RAREFACTION ANALYSIS ===")
    print(f"Total samples: {len(df)}")
    print(f"Min reads per sample: {df['frequency'].min():,}")
    print(f"Max reads per sample: {df['frequency'].max():,}")
    print(f"Mean reads per sample: {df['frequency'].mean():.0f}")
    print(f"Median reads per sample: {df['frequency'].median():.0f}")
    
    # Calculate percentiles
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    print("\nPercentiles:")
    for p in percentiles:
        value = df['frequency'].quantile(p/100)
        print(f"  {p}th percentile: {value:,.0f}")
    
    # Suggest rarefaction depths
    print("\nSuggested rarefaction depths:")
    conservative = df['frequency'].quantile(0.10)  # 10th percentile
    moderate = df['frequency'].quantile(0.25)      # 25th percentile
    aggressive = df['frequency'].quantile(0.50)    # 50th percentile
    
    print(f"  Conservative (retain 90% samples): {conservative:,.0f}")
    print(f"  Moderate (retain 75% samples): {moderate:,.0f}")
    print(f"  Aggressive (retain 50% samples): {aggressive:,.0f}")
    
    # Count samples at different thresholds
    thresholds = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 40000, 50000]
    print("\nSamples retained at different thresholds:")
    for threshold in thresholds:
        retained = (df['frequency'] >= threshold).sum()
        percent = (retained / len(df)) * 100
        print(f"  {threshold:,} reads: {retained} samples ({percent:.1f}%)")

if __name__ == "__main__":
    # Analyze 16S data
    qzv_16s = "qiime2/filtered/CFM_16S_final_filtered_summary.qzv"
    if os.path.exists(qzv_16s):
        df_16s = extract_sample_stats(qzv_16s)
        analyze_rarefaction_depth(df_16s, "16S BACTERIAL")
    else:
        print(f"16S summary file not found: {qzv_16s}")
    
    # Analyze ITS1 data
    qzv_its1 = "qiime2/filtered/CFM_ITS1_final_filtered_summary.qzv"
    if os.path.exists(qzv_its1):
        df_its1 = extract_sample_stats(qzv_its1)
        analyze_rarefaction_depth(df_its1, "ITS1 FUNGAL")
    else:
        print(f"ITS1 summary file not found: {qzv_its1}")
