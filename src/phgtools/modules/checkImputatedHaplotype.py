#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module to check the percentage of ranges in an imputed hVCF file that match 
source pangenome hVCF files from a PHGv2 database.

Reads hVCF files directly.
Optimized using set operations for fast key matching.
"""

import os
import glob
import gzip
import matplotlib.pyplot as plt
import argparse
import sys
from tqdm import tqdm
import pandas as pd
from datetime import datetime


def list_pangenome_hvcfs(hvcf_folder):
    """
    Finds and returns a list of gzipped hVCF files in the specified folder.
    """
    hvcf_files = glob.glob(os.path.join(hvcf_folder, '*.h.vcf.gz'))
    if not hvcf_files:
        print(f"Error: No *.h.vcf.gz files found in the folder: {hvcf_folder}")
        sys.exit(1)
        
    print(f"There are {len(hvcf_files)} pangenome source genomes to check.")
    return hvcf_files


def parse_hvcf_keys(hvcf_path, verbose=False):
    """
    Parse hVCF file and extract unique haplotype keys efficiently.
    
    Args:
        hvcf_path (str): Path to hVCF file (gzipped or plain text)
        verbose (bool): If True, print processing messages
        
    Returns:
        set: Set of unique haplotype keys
    """
    if hvcf_path.endswith('.gz'):
        open_func = gzip.open(hvcf_path, 'rt')
    else:
        open_func = open(hvcf_path, 'r')
    
    keys_set = set()
    try:
        with open_func as f:
            for line in f:
                if line.startswith("#"):
                    continue
                
                # VCF fields: CHROM POS ID REF ALT ...
                fields = line.strip().split("\t")
                if len(fields) > 4:
                    # The key (haplotype ID) is in the ALT column (index 4)
                    key_field = fields[4].strip()
                    
                    # PHG hVCF uses comma-separated keys in ALT, e.g., <key1>,<key2>
                    for k in key_field.split(','):
                        k = k.strip('<>').strip()
                        if k and k != ".":
                            keys_set.add(k)
    except Exception as e:
        if verbose:
            print(f"Error reading file {hvcf_path}: {e}")
        return set()
    
    return keys_set


def parse_hvcf_with_ranges(hvcf_path, verbose=False):
    """
    Parse hVCF file and extract haplotype keys with their range information.
    Calculates base pairs covered for each haplotype key from ##ALT header lines.
    
    Args:
        hvcf_path (str): Path to hVCF file (gzipped or plain text)
        verbose (bool): If True, print processing messages
        
    Returns:
        dict: Dictionary mapping haplotype keys to their span (end - start) in base pairs
    """
    if hvcf_path.endswith('.gz'):
        open_func = gzip.open(hvcf_path, 'rt')
    else:
        open_func = open(hvcf_path, 'r')
    
    keys_bp = {}
    try:
        with open_func as f:
            for line in f:
                # Parse header lines containing ALT definitions with range info
                if line.startswith("##ALT="):
                    # Extract the haplotype ID and RefRange from the ##ALT line
                    # Format: ##ALT=<ID=checksum,...,RefRange=chr:start-end,...>
                    import re
                    
                    # Extract ID (haplotype key)
                    id_match = re.search(r'ID=([^,>]+)', line)
                    if id_match:
                        hap_id = id_match.group(1).strip()
                        
                        # Extract RefRange
                        refrange_match = re.search(r'RefRange=([^:]+):(\d+)-(\d+)', line)
                        if refrange_match:
                            start = int(refrange_match.group(2))
                            end = int(refrange_match.group(3))
                            range_bp = max(0, end - start)
                            keys_bp[hap_id] = range_bp
    except Exception as e:
        if verbose:
            print(f"Error reading file {hvcf_path}: {e}")
        return {}
    
    return keys_bp




def check_imputation(source_hvcf_list, imputed_hvcf_path, verbose=False):
    """
    Optimized comparison of haplotype keys using set intersection.
    Compares the haplotype keys in the imputed hVCF against all source hVCFs
    to calculate the percentage of matching ranges and base pairs covered.
    
    Args:
        source_hvcf_list (list): Paths to the pangenome source hVCF files.
        imputed_hvcf_path (str): Path to the single imputed hVCF file.
        verbose (bool): If True, print detailed processing messages.

    Returns:
        pd.DataFrame: DataFrame with columns ['Genome', 'Matches', 'Total', 'Match_Percentage',
                     'Imputed_BP_Covered', 'Matched_BP_Covered', 'Matched_BP_Percentage']
    """

    # 1. Extract imputed hVCF keys into a set for O(1) lookup, and get base pair information
    if verbose:
        print(f"Processing imputed hVCF: {imputed_hvcf_path}")
    
    imputed_keys_set = parse_hvcf_keys(imputed_hvcf_path, verbose=verbose)
    imputed_keys_bp = parse_hvcf_with_ranges(imputed_hvcf_path, verbose=verbose)
    total_ranges = len(imputed_keys_set)
    total_imputed_bp = sum(imputed_keys_bp.values())
    
    if verbose:
        print(f"Total unique haplotype IDs in imputed file: {total_ranges}")
        print(f"Total base pairs covered in imputed file: {total_imputed_bp}")
    
    if total_ranges == 0:
        print("Error: No valid haplotype ranges found in the imputed file.")
        return pd.DataFrame()

    # 2. Compare against source hVCFs using set intersection (O(n) instead of O(n²))
    if verbose:
        print("\nStarting comparison against pangenome source genomes...")

    results = []
    
    for source_hvcf_path in tqdm(source_hvcf_list, desc="Processing source genomes", disable=not verbose):
        file_name = os.path.basename(source_hvcf_path).split(".")[0]
        
        try:
            # Extract keys from source file
            source_keys_set = parse_hvcf_keys(source_hvcf_path, verbose=verbose)
            source_keys_bp = parse_hvcf_with_ranges(source_hvcf_path, verbose=verbose)
            
            # Set intersection: find matching keys (vectorized, very fast)
            matched_keys = imputed_keys_set.intersection(source_keys_set)
            match_count = len(matched_keys)
            
            # Calculate base pairs covered by matched haplotypes
            matched_bp = sum(imputed_keys_bp.get(key, 0) for key in matched_keys)
            
            match_percentage = round(match_count / total_ranges * 100, 2)
            matched_bp_percentage = round(matched_bp / total_imputed_bp * 100, 2) if total_imputed_bp > 0 else 0.0
            
            if verbose:
                print(f"  {file_name}: {match_count}/{total_ranges} matches ({match_percentage}%), "
                      f"{matched_bp}/{total_imputed_bp} BP covered ({matched_bp_percentage}%)")
            
            results.append({
                'Genome': file_name,
                'Matches': match_count,
                'Total': total_ranges,
                'Match_Percentage': match_percentage,
                'Imputed_BP_Covered': total_imputed_bp,
                'Matched_BP_Covered': matched_bp,
                'Matched_BP_Percentage': matched_bp_percentage
            })
            
        except Exception as e:
            if verbose:
                print(f"  Error processing {file_name}: {e}")
            continue

    # Convert results to DataFrame for better handling and export
    results_df = pd.DataFrame(results)
    return results_df


def plot_imputation(results_df, output_dir=".", verbose=False):
    """
    Generates a bar plot of the imputation match percentages and BP coverage, and saves it.
    
    Args:
        results_df (pd.DataFrame): DataFrame with match results.
        output_dir (str): Directory to save output files.
        verbose (bool): If True, print file locations.
    """
    if results_df.empty:
        print("No data available to plot.")
        return None

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Sort by match percentage in descending order
    results_df = results_df.sort_values('Match_Percentage', ascending=False)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(10, len(results_df) * 0.5), 10), sharex=True)
    
    # Plot 1: Haplotype Match Percentage
    ax1.bar(results_df['Genome'], results_df['Match_Percentage'], color='skyblue', edgecolor='navy', alpha=0.7)
    ax1.set_ylabel('Haplotype Match (%)', fontsize=12, fontweight='bold')
    ax1.set_title('Imputed Haplotype Match Percentage per Source Genome', fontsize=14, fontweight='bold')
    ax1.set_ylim(0, 105)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Plot 2: Base Pair Coverage Percentage
    ax2.bar(results_df['Genome'], results_df['Matched_BP_Percentage'], color='lightgreen', edgecolor='darkgreen', alpha=0.7)
    ax2.set_ylabel('Base Pair Match (%)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Source Genome', fontsize=12, fontweight='bold')
    ax2.set_title('Imputed Base Pair Coverage Percentage per Source Genome', fontsize=14, fontweight='bold')
    ax2.set_ylim(0, 105)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.tight_layout()
    
    # Save figure
    plot_file = os.path.join(output_dir, 'imputation_match_percentage.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    if verbose:
        print(f"Plot saved to: {plot_file}")
    
    return plot_file


def main(args=None):
    """
    Main function to handle argument parsing and execution flow.
    
    Args:
        args: List of command line arguments (for testing/module integration)
    """
    
    # Use provided args or sys.argv
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="check-imputated-haplotype",
        description="Check the percentage of ranges in an imputed hVCF that match source pangenome hVCFs."
    )
    
    parser.add_argument("pangenome_folder", help="Path to the folder containing source pangenome *.h.vcf.gz files.")
    parser.add_argument("imputed_hvcf", help="Path to the imputed hVCF file (can be gzipped or not).")
    parser.add_argument("-o", "--output", type=str, default=".", help="Output directory for results (default: current directory)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print detailed processing messages.")
    
    parsed_args = parser.parse_args(args)

    # Get the list of source hVCFs
    source_hvcf_list = list_pangenome_hvcfs(parsed_args.pangenome_folder)

    # Perform the check and get results as DataFrame
    results_df = check_imputation(source_hvcf_list, parsed_args.imputed_hvcf, verbose=parsed_args.verbose)

    if results_df.empty:
        print("No results to display.")
        return

    # Plot the results
    plot_file = plot_imputation(results_df, output_dir=parsed_args.output, verbose=parsed_args.verbose)

    # Save results to TSV
    tsv_file = os.path.join(parsed_args.output, 'imputation_match_results.tsv')
    results_df.to_csv(tsv_file, sep='\t', index=False)

    # Print summary
    print("\n" + "="*90)
    print("IMPUTATION MATCH SUMMARY")
    print("="*90)
    print(f"\nTotal genomes analyzed: {len(results_df)}")
    print(f"\nMatch Percentage Statistics:")
    print(f"  Mean match percentage: {results_df['Match_Percentage'].mean():.2f}%")
    print(f"  Min match percentage:  {results_df['Match_Percentage'].min():.2f}%")
    print(f"  Max match percentage:  {results_df['Match_Percentage'].max():.2f}%")
    print(f"\nBase Pairs Covered Statistics:")
    print(f"  Mean matched BP percentage: {results_df['Matched_BP_Percentage'].mean():.2f}%")
    print(f"  Min matched BP percentage:  {results_df['Matched_BP_Percentage'].min():.2f}%")
    print(f"  Max matched BP percentage:  {results_df['Matched_BP_Percentage'].max():.2f}%")
    print(f"\nDetailed Results:")
    print(results_df.to_string(index=False))
    print("\n" + "="*90)
    print("OUTPUT FILES")
    print("="*90)
    print(f"TSV Results:  {tsv_file}")
    if plot_file:
        print(f"Plot Image:   {plot_file}")
    print("="*90 + "\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nScript interrupted by user. Exiting...")
        sys.exit(0)
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        sys.exit(1)
