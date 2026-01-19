#!/usr/bin/env python
# coding: utf-8

"""
Detect core, unique, and accessory ranges from a pangenome hVCF file.

This module analyzes hVCF files from PHGv2 to identify:
- Core ranges: Present in all genomes
- Unique ranges: Present in only one genome
- Accessory ranges: Present in 2 to (n-1) genomes

Where n is the total number of genomes in the pangenome.
"""

import os
import argparse
import sys
from pathlib import Path
import gzip
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def IsValidHVCFFile(file_path):
    """
    Check if a file has a valid hVCF extension.
    
    Parameters:
    -----------
    file_path : str
        Path to file to check
        
    Returns:
    --------
    bool : True if file has valid hVCF extension
    """
    return str(file_path).endswith(('.h.vcf', '.h.vcf.gz'))


def NumberOfGenomes(hvcf_input):
    """
    Determine the maximum number of genomes in the pangenome from hVCF file.
    
    Parameters:
    -----------
    hvcf_input : str
        Path to the merged pangenome hVCF file
    
    Returns:
    --------
    int : Maximum number of keys in any line (maximum number of genomes)
    """
    merged_file = hvcf_input
    
    num_of_genomes = 0
    
    # Determine if file is gzipped
    is_gzipped = str(merged_file).endswith('.gz')
    
    if is_gzipped:
        file_handle = gzip.open(merged_file, 'rt')
    else:
        file_handle = open(merged_file, 'r')
    
    try:
        for line in file_handle:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) > 4:
                keys = fields[4].split(',')
                if len(keys) > num_of_genomes:
                    num_of_genomes = len(keys)
    finally:
        file_handle.close()
    
    return num_of_genomes


def DetermineCoreRanges(hvcf_input):
    """
    Analyze pangenome hVCF file to determine core, unique, and accessory ranges.
    
    Parameters:
    -----------
    hvcf_input : str
        Path to the merged pangenome hVCF file
    
    Returns:
    --------
    tuple : (total_ranges, unique_ranges, core_ranges, accessory_count, dict_ranges_count, genotype_counts)
        - total_ranges: Total number of ranges in pangenome
        - unique_ranges: Ranges present in only 1 genome
        - core_ranges: Ranges present in all genomes
        - accessory_count: Ranges present in 2 to (n-1) genomes
        - dict_ranges_count: Dictionary with count by genome number
        - genotype_counts: Dictionary counting genotypes in ranges with exactly 2 genomes
    """
    merged_file = hvcf_input
    
    num_of_genomes = NumberOfGenomes(merged_file)
    genome_numbers = list(range(1, num_of_genomes + 1))
    
    # Determine if file is gzipped
    is_gzipped = str(merged_file).endswith('.gz')
    
    # First pass: count totals
    total_ranges = 0
    core_ranges = 0
    unique_ranges = 0
    genotype_counts = {}
    
    if is_gzipped:
        file_handle = gzip.open(merged_file, 'rt')
    else:
        file_handle = open(merged_file, 'r')
    
    try:
        for line in file_handle:
            if line.startswith('#'):
                continue
            total_ranges += 1
            fields = line.strip().split('\t')
            if len(fields) > 4:
                keys = fields[4].split(',')
                if len(keys) == num_of_genomes:
                    core_ranges += 1
                if len(keys) == 1:
                    unique_ranges += 1
                # Track genotypes for ranges with exactly 2 genomes
                if len(keys) == 2:
                    for key in keys:
                        genotype_counts[key] = genotype_counts.get(key, 0) + 1
    finally:
        file_handle.close()
    
    # Second pass: count by genome number
    dict_ranges_count = {}
    for genome_number in genome_numbers:
        count = 0
        
        if is_gzipped:
            file_handle = gzip.open(merged_file, 'rt')
        else:
            file_handle = open(merged_file, 'r')
        
        try:
            for line in file_handle:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) > 4:
                    keys = fields[4].split(',')
                    if len(keys) == genome_number:
                        count += 1
        finally:
            file_handle.close()
        
        dict_ranges_count[genome_number] = count
    
    # Calculate accessory
    accessory_count = 0
    for genome_number in genome_numbers:
        if genome_number != 1 and genome_number != num_of_genomes:
            accessory_count += dict_ranges_count[genome_number]
    
    # Validate
    if total_ranges != core_ranges + unique_ranges + accessory_count:
        # Debug print statements
        print("\n" + "="*70)
        print("ERROR: Validation failed")
        print("="*70)
        print(f"Total ranges: {total_ranges}")
        print(f"Core ranges: {core_ranges}")
        print(f"Unique ranges: {unique_ranges}")
        print(f"Accessory ranges: {accessory_count}")
        print(f"Sum of categories: {core_ranges + unique_ranges + accessory_count}")
        print(f"\nNumber of genomes: {num_of_genomes}")
        print(f"Ranges by genome count:")
        for gn in sorted(dict_ranges_count.keys()):
            print(f"  {gn}: {dict_ranges_count[gn]}")
        print("="*70 + "\n")
        raise ValueError('The sum of core, unique and accessory ranges does not match the total')
    
    return total_ranges, unique_ranges, core_ranges, accessory_count, dict_ranges_count, genotype_counts


def PrintCoreRangesStats(hvcf_input):
    """
    Print formatted statistics about core, unique, and accessory ranges,
    and generate visualizations (bar plot and pie chart).
    
    Parameters:
    -----------
    hvcf_input : str
        Path to the merged pangenome hVCF file
    """
    merged_file = hvcf_input
    
    total_ranges, unique_ranges, core_ranges, accessory_count, dict_ranges_count, genotype_counts = DetermineCoreRanges(merged_file)
    num_of_genomes = NumberOfGenomes(merged_file)
    
    print("\n" + "="*70)
    print("PANGENOME RANGE ANALYSIS")
    print("="*70)
    print(f"Input file: {hvcf_input}")
    print(f"Total number of ranges: {total_ranges}")
    print(f"Number of genomes: {num_of_genomes}")
    print(f"Core ranges (present in all genomes): {core_ranges} ({core_ranges/total_ranges*100:.2f}%)")
    print(f"Unique ranges (present in 1 genome): {unique_ranges} ({unique_ranges/total_ranges*100:.2f}%)")
    print(f"Accessory ranges (present in 2-{num_of_genomes-1} genomes): {accessory_count} ({accessory_count/total_ranges*100:.2f}%)")
    print("="*70)
    print("\nRanges by genome count:")
    print("-" * 70)
    for genome_number in sorted(dict_ranges_count.keys()):
        count = dict_ranges_count[genome_number]
        percentage = (count / total_ranges * 100) if total_ranges > 0 else 0
        range_type = ""
        if genome_number == 1:
            range_type = " (Unique)"
        elif genome_number == num_of_genomes:
            range_type = " (Core)"
        else:
            range_type = " (Accessory)"
        print(f"Ranges present in {genome_number:>2d} genome(s): {count:>8d} ({percentage:>6.2f}%){range_type}")
    print("-" * 70 + "\n")
    
    # Generate visualizations
    GenerateVisualizations(dict_ranges_count, core_ranges, unique_ranges, accessory_count, 
                          total_ranges, num_of_genomes, hvcf_input)


def GenerateVisualizations(dict_ranges_count, core_ranges, unique_ranges, accessory_count, 
                          total_ranges, num_of_genomes, hvcf_input):
    """
    Generate bar plot and pie chart visualizations for range analysis.
    
    Parameters:
    -----------
    dict_ranges_count : dict
        Dictionary with count by genome number
    core_ranges : int
        Number of core ranges
    unique_ranges : int
        Number of unique ranges
    accessory_count : int
        Number of accessory ranges
    total_ranges : int
        Total number of ranges
    num_of_genomes : int
        Total number of genomes
    hvcf_input : str
        Path to input file (used for title)
    """
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Bar plot: Ranges by genome count
    genome_numbers = sorted(dict_ranges_count.keys())
    counts = [dict_ranges_count[gn] for gn in genome_numbers]
    colors = ['#E74C3C' if gn == 1 else '#3498DB' if gn == num_of_genomes else '#F39C12' 
              for gn in genome_numbers]
    
    ax1.bar(genome_numbers, counts, color=colors, edgecolor='black', linewidth=1.2)
    ax1.set_xlabel('Number of Genomes', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Ranges', fontsize=12, fontweight='bold')
    ax1.set_title('Ranges by Genome Count', fontsize=14, fontweight='bold')
    ax1.set_xticks(genome_numbers)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add legend for bar plot
    unique_patch = mpatches.Patch(facecolor='#E74C3C', label='Unique (1 genome)', edgecolor='black', linewidth=1.2)
    core_patch = mpatches.Patch(facecolor='#3498DB', label=f'Core ({num_of_genomes} genomes)', edgecolor='black', linewidth=1.2)
    accessory_patch = mpatches.Patch(facecolor='#F39C12', label='Accessory (2-12 genomes)', edgecolor='black', linewidth=1.2)
    ax1.legend(handles=[unique_patch, core_patch, accessory_patch], loc='upper left')
    
    # Pie chart: Core vs Unique vs Accessory
    labels = [f'Core\n{core_ranges:,}\n({core_ranges/total_ranges*100:.1f}%)',
              f'Accessory\n{accessory_count:,}\n({accessory_count/total_ranges*100:.1f}%)',
              f'Unique\n{unique_ranges:,}\n({unique_ranges/total_ranges*100:.1f}%)']
    sizes = [core_ranges, accessory_count, unique_ranges]
    colors_pie = ['#3498DB', '#F39C12', '#E74C3C']
    explode = (0.05, 0.02, 0.05)
    
    ax2.pie(sizes, labels=labels, colors=colors_pie, autopct='', startangle=90,
            explode=explode, textprops={'fontsize': 11, 'fontweight': 'bold'},
            wedgeprops=dict(edgecolor='black', linewidth=1.5))
    
    # Save figure
    output_file = hvcf_input.replace('.h.vcf.gz', '').replace('.h.vcf', '') + '_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to: {output_file}\n")
    
    # Display plot
    plt.show()


def main(args):
    """
    Main function to analyze pangenome core ranges.
    
    Parameters:
    -----------
    args : list
        Command line arguments
    """
    parser = argparse.ArgumentParser(
        prog="phgtools core-range-detector",
        description="Analyze core, unique, and accessory ranges from merged pangenome hVCF file."
    )
    
    parser.add_argument(
        "merged_hvcf",
        metavar="HVCF_FILE",
        help="Path to merged pangenome hVCF file (.h.vcf or .h.vcf.gz)"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    
    parsed_args = parser.parse_args(args)
    
    hvcf_input = parsed_args.merged_hvcf
    verbose = parsed_args.verbose
    
    try:
        if verbose:
            # print(f"Analyzing merged pangenome file: {hvcf_input}\n")
            pass
        
        # Print statistics
        PrintCoreRangesStats(hvcf_input)
    
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
