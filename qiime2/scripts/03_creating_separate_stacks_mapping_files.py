#!/usr/bin/env python3
"""
Script to generate stacks combinatorial tag mapping files from TSV file originating from metadata
I run it using python/3.12.6 module
Input: TSV file with sample metadata and tag sequences: stacks_sample_mapping_all_sublibraries.txt
Output: Individual TSV files for each sublibrary (48 files)
"""

import os

def create_mapping_files(input_file, output_dir):
    """
    Read TSV file and create individual mapping files for each sublibrary
    
    Parameters:
    input_file (str): Path to TSV file with my data
    output_dir (str): Directory to save the mapping files
    """
    
    print(f" Reading file: {input_file}")
    
    # Read the TSV file
    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f" File not found: {input_file}")
        return
    
    print(f" Read {len(lines)} lines from file")
    
    # Remove header line and empty lines
    header = lines[0].strip().split('\t')
    data_lines = [line.strip() for line in lines[1:] if line.strip()]
    
    print(f" Header columns: {header}")
    print(f" Data rows: {len(data_lines)}")
    
    # Column positions based on TSV file:
    # 0=sample-id, 1=sublibrary_id, 2=16s_forward_primer, 3=16s_forward_primer_tag, 
    # 4=16s_reverse_primer, 5=16s_reverse_primer_tag, 6=its1_forward_primer, 
    # 7=its1_forward_primer_tag, 8=its1_reverse_primer, 9=its1_reverse_primer_tag
    SAMPLE_ID = 0           # sample-id  
    SUBLIBRARY_ID = 1       # sublibrary_id
    S16_FORWARD_TAG = 3     # 16s_forward_primer_tag
    S16_REVERSE_TAG = 5     # 16s_reverse_primer_tag  
    ITS1_FORWARD_TAG = 7    # its1_forward_primer_tag
    ITS1_REVERSE_TAG = 9    # its1_reverse_primer_tag
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f" Output directory: {output_dir}")
    
    # Group data by sublibrary
    sublibraries = {}
    
    for line in data_lines:
        columns = line.split('\t')
        
        # Skip lines that don't have enough columns
        if len(columns) < 10:
            print(f" Skipping line (not enough columns): {line[:50]}...")
            continue
        
        sample_id = columns[SAMPLE_ID]
        sublibrary_id = columns[SUBLIBRARY_ID]
        
        # Skip empty rows
        if not sample_id or not sublibrary_id:
            continue
        
        # Add to sublibrary group
        if sublibrary_id not in sublibraries:
            sublibraries[sublibrary_id] = []
        
        # Get tag sequences
        s16_forward = columns[S16_FORWARD_TAG]
        s16_reverse = columns[S16_REVERSE_TAG]
        its1_forward = columns[ITS1_FORWARD_TAG]
        its1_reverse = columns[ITS1_REVERSE_TAG]
        
        # Add both 16S and ITS1 entries for this sample
        sublibraries[sublibrary_id].append({
            'forward_tag': s16_forward,
            'reverse_tag': s16_reverse,
            'sample_name': f"{sample_id}_16S"
        })
        
        sublibraries[sublibrary_id].append({
            'forward_tag': its1_forward,
            'reverse_tag': its1_reverse,
            'sample_name': f"{sample_id}_ITS1"
        })
    
    print(f" Found {len(sublibraries)} unique sublibraries")
    
    # Create mapping file for each sublibrary
    files_created = 0
    
    for sublibrary_id in sorted(sublibraries.keys()):  # Sort for consistent output
        samples = sublibraries[sublibrary_id]
        
        print(f"\n Processing sublibrary: {sublibrary_id}")
        print(f"   Samples: {len(samples)//2} ({len(samples)} tag combinations)")
        
        # Create output file
        output_file = os.path.join(output_dir, f"{sublibrary_id}_tags.tsv")
        
        with open(output_file, 'w') as file:
            for sample in samples:
                # Write in Stacks format: forward_tag\treverse_tag\tsample_name
                file.write(f"{sample['forward_tag']}\t{sample['reverse_tag']}\t{sample['sample_name']}\n")
        
        print(f" Created: {output_file}")
        
        # Show first few lines as example
        print(f" Example entries:")
        for i, sample in enumerate(samples[:4]):  # Show first 4 entries
            print(f"      {sample['forward_tag']}\t{sample['reverse_tag']}\t{sample['sample_name']}")
        
        files_created += 1
    
    print(f"\n Process complete!")
    print(f" Created {files_created} mapping files")
    print(f" Ready for Stacks process_radtags!")

def main():
    """
    Main function - update these paths for your setup
    """
    
    # INPUT: Path to your TSV file (UPDATE THIS)
    input_file = "qiime2/import/demux/stacks_sample_mapping_all_sublibraries.txt"  # Actual file is a subset of cacao_microbiome_metadata.xlsx with primer names substituted for tag sequences using =VLOOKUP() formula in Excel
    
    # OUTPUT: Directory for mapping files  
    output_dir = "qiime2/import/demux/internal_tag_mappings"
    
    # Create the mapping files
    create_mapping_files(input_file, output_dir)

if __name__ == "__main__":
    main()
