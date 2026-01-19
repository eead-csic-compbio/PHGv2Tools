#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to analyze pangenome evolution by calculating the cumulative growth of 
unique ranges (pangenome size) and unique allele keys (variation size) 
directly from a pre-computed hapIDranges.tsv file.
"""

import os
import glob
import sys
import subprocess
import matplotlib.pyplot as plt
import argparse
import gzip
from matplotlib.ticker import MaxNLocator

# --- Core Analysis Function ---

def analyze_pangenome_growth(hapid_ranges_file, reference_genome=None, verbose=False):
    """
    Analyzes the cumulative increase in unique haplotype keys and covered ranges 
    by sequentially "adding" samples from the TSV file.

    Args:
        hapid_ranges_file (str): Path to the hapIDranges.tsv file.
        reference_genome (str): Name of the reference genome. If None, user will be prompted.
        verbose (bool): If True, print detailed execution logs.

    Returns:
        tuple: (cumulative_ranges_dict, cumulative_keys_dict, samples)
    """
    
    # Data structures to store results
    # cumulative_ranges: {genome_count: total_ranges_with_at_least_one_haplotype_id}
    # cumulative_keys:   {genome_count: total_unique_haplotype_keys_found_so_far}
    cumulative_ranges = {} 
    cumulative_keys = {}   
    
    # Set to track all unique keys found across all samples added so far
    unique_keys_set = set()
    
    # Set to track the index of every range (row) that has at least one haplotype ID
    covered_ranges_indices = set() 
    
    range_data = [] # List of lists: [[Sample1_Key_R1, Sample2_Key_R1, ...], [Sample1_Key_R2, ...], ...]
    samples = []
    
    # 1. Read and Parse the TSV file
    if hapid_ranges_file.endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt')
    else:
        open_func = lambda f: open(f, 'r')
        
    print(f"Reading and parsing {hapid_ranges_file}...")

    try:
        with open_func(hapid_ranges_file) as f:
            # Parse Header
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().split('\t')
                    # Samples start from index 3 (after #CHROM, START, END)
                    samples = header[3:]
                    num_samples = len(samples)
                    break
            else:
                raise Exception("TSV file is empty or missing the #CHROM header line.")
                
            # Parse Range Data
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                hap_ids = fields[3:] 
                
                if len(hap_ids) != num_samples:
                     if verbose:
                         print(f"Skipping line {fields[0]}:{fields[1]}-{fields[2]} due to mismatch in column count.")
                     continue
                     
                range_data.append(hap_ids)
                
    except Exception as e:
        print(f"Error during file reading: {e}")
        sys.exit(1)

    # 2. Handle Reference Genome Selection
    if reference_genome is None:
        print("\n--- Available Genomes ---")
        for i, genome in enumerate(samples, 1):
            print(f"  {i}. {genome}")
        print()
        
        while True:
            try:
                user_input = input(f"Enter the reference genome name (or number 1-{num_samples}): ").strip()
                
                # Check if input is a number
                if user_input.isdigit():
                    idx = int(user_input) - 1
                    if 0 <= idx < num_samples:
                        reference_genome = samples[idx]
                        break
                    else:
                        print(f"Error: Number must be between 1 and {num_samples}. Try again.")
                else:
                    # Check if input is a genome name
                    if user_input in samples:
                        reference_genome = user_input
                        break
                    else:
                        print(f"Error: '{user_input}' not found in genome list. Try again.")
            except KeyboardInterrupt:
                print("\nAborted by user.")
                sys.exit(0)
    
    # 3. Validate reference genome
    if reference_genome not in samples:
        print(f"Error: Reference genome '{reference_genome}' not found in the file.")
        sys.exit(1)
    
    if verbose:
        print(f"\nReference genome: {reference_genome}\n")

    # 4. Count ranges per genome and sort for ranges amplification analysis
    ranges_per_genome = {}
    for sample in samples:
        ranges_per_genome[sample] = 0
    
    for hap_ids in range_data:
        for sample_index, key in enumerate(hap_ids):
            key_clean = key.strip('<>').strip()
            if key_clean and key_clean != ".":
                ranges_per_genome[samples[sample_index]] += 1
    
    # Sort genomes by number of ranges (ascending), with reference genome at the end
    sorted_samples_by_ranges = sorted(
        [s for s in samples if s != reference_genome],
        key=lambda x: ranges_per_genome[x]
    )
    sorted_samples_by_ranges.append(reference_genome)
    
    # Reorder range_data columns to match sorted_samples_by_ranges
    sample_to_index = {sample: i for i, sample in enumerate(samples)}
    sorted_indices = [sample_to_index[sample] for sample in sorted_samples_by_ranges]
    
    reordered_range_data = []
    for hap_ids in range_data:
        reordered_range_data.append([hap_ids[i] for i in sorted_indices])
    
    range_data = reordered_range_data
    samples = sorted_samples_by_ranges
    
    if verbose:
        print(f"Genomes sorted by ranges (ascending to reference):")
        for sample in samples:
            print(f"  {sample}: {ranges_per_genome[sample]} ranges")
    
    # Check if reference has the highest number of ranges
    reference_range_count = ranges_per_genome[reference_genome]
    max_range_count = max(ranges_per_genome.values())
    
    if reference_range_count < max_range_count:
        print("\n" + "!"*90)
        print("WARNING: The reference genome does NOT have the highest number of ranges!")
        print(f"  Reference: {reference_genome} has {reference_range_count} ranges")
        print(f"  Maximum:   {[g for g, c in ranges_per_genome.items() if c == max_range_count][0]} has {max_range_count} ranges")
        print("  Are you sure you indicated the proper pangenome reference?")
        print("!"*90 + "\n")
    
    # 5. Analyze Cumulative Growth
    if not range_data:
        print("Error: No valid range data found in the file.")
        return {}, {}, []
        
    print(f"Analysis based on {len(range_data)} total ranges and {num_samples} samples.")
    
    # Track new vs existing keys for each sample
    keys_per_sample = []
    ranges_per_sample = []
    
    # Iterate column by column (sample by sample)
    for sample_index in range(num_samples):
        sample_name = samples[sample_index]
        new_keys_added = 0
        shared_keys = 0
        keys_in_this_sample = set()
        ranges_in_this_sample = 0
        
        # Iterate row by row (range by range)
        for range_index, hap_ids in enumerate(range_data):
            # Extract the key for the current sample in this range
            key = hap_ids[sample_index].strip('<>').strip()
            
            if key and key != ".":
                ranges_in_this_sample += 1
                
                # Only count each key once per sample
                if key not in keys_in_this_sample:
                    keys_in_this_sample.add(key)
                    
                    # Check if this is a new key to the pangenome
                    if key not in unique_keys_set:
                        new_keys_added += 1
                        unique_keys_set.add(key)
                    else:
                        # Key already exists from previous genomes
                        shared_keys += 1
                
                # Mark range as covered (Size)
                covered_ranges_indices.add(range_index)
        
        # Store cumulative results
        genome_count = sample_index + 1
        cumulative_keys[genome_count] = len(unique_keys_set)
        cumulative_ranges[genome_count] = len(covered_ranges_indices)
        
        keys_per_sample.append({
            'sample': sample_name,
            'new_keys': new_keys_added,
            'shared_keys': shared_keys,
            'total_cumulative': len(unique_keys_set)
        })
        
        ranges_per_sample.append({
            'sample': sample_name,
            'ranges_in_genome': ranges_in_this_sample,
            'cumulative_ranges': len(covered_ranges_indices)
        })
        
        if verbose:
            print(f"\n--- Sample: {sample_name} ({genome_count}/{num_samples}) ---")
            print(f"  Ranges in this genome: {ranges_in_this_sample}")
            print(f"  Total Ranges Covered: {len(covered_ranges_indices)}")
            print(f"  Total Unique Keys:    {len(unique_keys_set)}")
            print(f"  New keys in this genome: {new_keys_added}")
            print(f"  Keys shared with previous genomes: {shared_keys}")

    # Print ranges amplification summary
    print("\n" + "="*90)
    print("SUMMARY: Ranges Amplification (Pangenome Size)")
    print("="*90)
    print(f"{'Genome':<20} {'Ranges in Genome':<20} {'Cumulative Ranges':<20}")
    print("-"*90)
    for item in ranges_per_sample:
        print(f"{item['sample']:<20} {item['ranges_in_genome']:<20} {item['cumulative_ranges']:<20}")
    print("="*90)

    # Print keys summary table
    print("\n" + "="*90)
    print("SUMMARY: Keys Variation (Pangenome Variation)")
    print("="*90)
    print(f"{'Genome':<20} {'New Keys':<15} {'Shared Keys':<18} {'Cumulative Total':<18}")
    print("-"*90)
    for item in keys_per_sample:
        print(f"{item['sample']:<20} {item['new_keys']:<15} {item['shared_keys']:<18} {item['total_cumulative']:<18}")
    print("="*90 + "\n")

    return cumulative_ranges, cumulative_keys, samples


# --- Plotting Function ---

def PlotResults(hapid_ranges_file, ranges_dict, keys_dict, samples):
    """
    Plots the Ranges Amplification and Alleles Variation slopes with sample labels.
    
    Args:
        hapid_ranges_file (str): Path to the input TSV file.
        ranges_dict (dict): Cumulative ranges data.
        keys_dict (dict): Cumulative keys data.
        samples (list): List of sample names in order.
    """
    
    if not ranges_dict or not keys_dict:
        print("No data available to plot.")
        return
        
    # Determine the output folder
    output_folder = os.path.dirname(hapid_ranges_file) or '.'
    
    # --- Plot 1: Ranges Amplification Slope (Pangenome Size) ---
    int_keys = list(ranges_dict.keys())
    ranges_values = list(ranges_dict.values())
    
    plt.figure(figsize=(14, 8))
    plt.plot(int_keys, ranges_values, marker='o', linestyle='-', color='blue')

    # Add genome name labels with alternating positions to avoid overlap
    for i, txt in enumerate(ranges_values):
        genome_name = samples[i] if i < len(samples) else f"Sample_{i+1}"
        
        # Alternate position: above and below the point
        if i % 2 == 0:
            xytext = (0, 20)
            va = 'bottom'
        else:
            xytext = (0, -30)
            va = 'top'
        
        plt.annotate(genome_name, 
                    (int_keys[i], ranges_values[i]), 
                    textcoords="offset points", 
                    xytext=xytext, 
                    ha='center',
                    va=va,
                    fontsize=10,
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='yellow', alpha=0.7, edgecolor='gray'))

    plt.xlabel('Number of Genomes Sequentially Added', fontsize=12)
    plt.ylabel('Total Number of Ranges Covered', fontsize=12)
    plt.title('Ranges Amplification Slope (Pangenome Size Growth)', fontsize=14)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid(True, linestyle='--', alpha=0.6)
    
    range_plot_path = os.path.join(output_folder, "RangesAmplificationSlope_from_TSV.png")
    plt.savefig(range_plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"\nPlot saved: {range_plot_path}")

    # --- Plot 2: Alleles Variation Slope (Pangenome Variation) ---
    int_keys_var = list(keys_dict.keys())
    keys_values = list(keys_dict.values())

    plt.figure(figsize=(14, 8))
    plt.plot(int_keys_var, keys_values, marker='o', linestyle='-', color='red')
    
    # Add genome name labels with alternating positions to avoid overlap
    for i, txt in enumerate(keys_values):
        genome_name = samples[i] if i < len(samples) else f"Sample_{i+1}"
        
        # Alternate position: above and below the point
        if i % 2 == 0:
            xytext = (0, 20)
            va = 'bottom'
        else:
            xytext = (0, -30)
            va = 'top'
        
        plt.annotate(genome_name, 
                    (int_keys_var[i], keys_values[i]), 
                    textcoords="offset points", 
                    xytext=xytext, 
                    ha='center',
                    va=va,
                    fontsize=10,
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='yellow', alpha=0.7, edgecolor='gray'))
    
    plt.xlabel('Number of Genomes Sequentially Added', fontsize=12)
    plt.ylabel('Total Number of Unique Haplotype Keys', fontsize=12)
    plt.title('Ranges Variation Slope (Pangenome Variation Growth)', fontsize=14)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid(True, linestyle='--', alpha=0.6)
    
    key_plot_path = os.path.join(output_folder, "RangesVariationSlope_from_TSV.png")
    plt.savefig(key_plot_path, bbox_inches='tight', dpi=300)
    plt.close()

    print(f"Plot saved: {key_plot_path}")

def main(args=None):
    """
    Handles command-line arguments and orchestrates the analysis workflow.
    """
    if args is None:
        args = sys.argv[1:]
    
    parser = argparse.ArgumentParser(
        description="Analyze Pangenome size and variation growth by calculating cumulative metrics from a hapIDranges.tsv file."
    )
    
    parser.add_argument("hapid_ranges_file", 
                        help="Path to the hapIDranges.tsv (or .tsv.gz) file.")
    
    parser.add_argument("-r", "--reference", dest="reference_genome",
                        help="Name of the reference genome. If not provided, you will be prompted to select one.")
    
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Print detailed execution and cumulative growth logs.")
    
    parsed_args = parser.parse_args(args)
    
    hapid_ranges_file = parsed_args.hapid_ranges_file
    reference_genome = parsed_args.reference_genome
    verbose = parsed_args.verbose

    if not os.path.exists(hapid_ranges_file):
        raise FileNotFoundError(f"Input file not found: {hapid_ranges_file}")

    if verbose:
        print("\n--- SCRIPT ARGUMENTS ---")
        print(f"HapID Ranges File:  {hapid_ranges_file}")
        print(f"Reference Genome:   {reference_genome if reference_genome else 'Will be prompted'}")
        print("------------------------\n")

    # --- Workflow ---
    
    # 1. Analyze Growth
    ranges_dict, keys_dict, samples = analyze_pangenome_growth(hapid_ranges_file, reference_genome, verbose)
    
    # 2. Plot Results
    PlotResults(hapid_ranges_file, ranges_dict, keys_dict, samples)
    
    print("\nAnalysis complete.")

if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"\nFATAL ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        sys.exit(1)
