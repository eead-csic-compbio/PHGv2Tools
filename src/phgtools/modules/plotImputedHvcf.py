#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module to plot imputed hVCF files based on pangenomes haplotypes.

Converts imputed hVCF files to matplotlib-based haplotype painting plots,
showing genome coverage across chromosomes with colored blocks for each range.
"""

import os
import gzip
import re
import argparse
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import warnings
import time
from tqdm import tqdm

warnings.filterwarnings('ignore', category=UserWarning)


def assign_colors(genomes, reference_name):
    """
    Assign a color to each genome.
    Reference genome is assigned as Black (#000000).
    Others use the same palette as haplopainting.
    
    Args:
        genomes (list): List of genome names
        reference_name (str): Name of the reference genome
        
    Returns:
        dict: Dictionary mapping genome names to hex color codes
    """
    # Color palette - same as haplopainting
    palette = [
        '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffd92f', 
        '#a65628', '#f781bf', '#999999', '#00ced1', '#000075', '#ffb300', 
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', 
        '#e5c494', '#b3b3b3'
    ]
    
    colors_dict = {}
    
    # Assign reference to black
    colors_dict[reference_name] = "#000000"
    
    # Assign other genomes from palette
    palette_idx = 0
    for genome in genomes:
        if genome != reference_name:
            colors_dict[genome] = palette[palette_idx % len(palette)]
            palette_idx += 1
    
    return colors_dict


def parse_hvcf_to_dataframe(imputed_hvcf, verbose=False):
    """
    Parse hVCF file and extract haplotype ranges into a pandas DataFrame.
    Parses ##ALT header lines to extract RefRange information and genome (SampleName).
    
    Args:
        imputed_hvcf (str): Path to imputed hVCF file
        verbose (bool): If True, print processing messages
        
    Returns:
        pd.DataFrame: DataFrame with columns [chr, ref_start, ref_end, genome]
    """
    data = []
    
    if imputed_hvcf.endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt')
    else:
        open_func = open
    
    with open_func(imputed_hvcf) as hvcf_file:
        for line in hvcf_file:
            if line.startswith('##ALT=<ID'):
                # Extract RefRange field
                m_range = re.search(r'RefRange=([^,>]+)', line)
                if not m_range:
                    continue
                
                range_val = m_range.group(1)  # e.g., chr5H:26022892-26058044
                
                if ':' not in range_val or '-' not in range_val:
                    continue
                
                try:
                    chr_part = range_val.split(':')[0].strip().replace('"', '').replace("'", "")
                    start = int(range_val.split(':')[1].split('-')[0])
                    end = int(range_val.split(':')[1].split('-')[1])
                except (ValueError, IndexError):
                    if verbose:
                        print(f"WARNING: Could not parse range: {range_val}")
                    continue
                
                if start >= end:
                    continue
                
                # Extract SampleName for genome identification (this is the pangenome genome it matches)
                m_sample = re.search(r'SampleName=([^,>]+)', line)
                genome = m_sample.group(1).strip().replace(" ", "") if m_sample else "unknown"
                
                data.append({
                    'chr': chr_part,
                    'ref_start': start,
                    'ref_end': end,
                    'genome': genome
                })
    
    df = pd.DataFrame(data)
    return df


def plot_haplotype_painting_with_data(dataframe, colors_dict, output_folder, imputed_hvcf, verbose=False):
    """
    Create a single haplotype painting figure with all chromosomes stacked vertically.
    Each range is colored according to which pangenome genome it matches.
    
    Args:
        dataframe (pd.DataFrame): DataFrame with [chr, ref_start, ref_end, genome]
        colors_dict (dict): Dictionary mapping genome names to hex colors
        output_folder (str): Output folder for PNG file
        imputed_hvcf (str): Name for output file
        verbose (bool): If True, print processing messages
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    if dataframe.empty:
        print("ERROR: No data to plot!")
        return
    
    # Add color column based on genome
    dataframe['color'] = dataframe['genome'].map(colors_dict)
    
    # Get unique chromosomes in order
    chromosomes = sorted(dataframe['chr'].unique())
    
    base_name = os.path.basename(imputed_hvcf).split('.')[0]
    
    # Create figure with subplots for each chromosome (no space between)
    n_chromosomes = len(chromosomes)
    fig_height = max(8, n_chromosomes * 0.6)  # Dynamic height
    fig, axes = plt.subplots(n_chromosomes, 1, figsize=(14, fig_height), sharex=False)
    
    # If only one chromosome, axes is not a list
    if n_chromosomes == 1:
        axes = [axes]
    
    # Track global min/max for consistent x-axis
    global_min = dataframe['ref_start'].min()
    global_max = dataframe['ref_end'].max()
    
    # Plot each chromosome
    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        chrom_data = dataframe[dataframe['chr'] == chrom]
        
        if chrom_data.empty:
            continue
        
        # Thin bar height
        track_height = 0.4
        
        # Draw background track (no visible background, transparent)
        ax.add_patch(Rectangle((global_min, -track_height/2), 
                               global_max - global_min, track_height,
                               facecolor='white', edgecolor='none', zorder=1))
        
        # Draw haplotype ranges colored by genome (no borders)
        for _, row in chrom_data.iterrows():
            ax.add_patch(Rectangle((row['ref_start'], -track_height/2),
                                  row['ref_end'] - row['ref_start'],
                                  track_height,
                                  facecolor=row['color'],
                                  edgecolor='none',
                                  zorder=2))
        
        # Style chromosome axis - minimal ticks
        ax.set_xlim(global_min, global_max)
        ax.set_ylim(-0.3, 0.3)
        ax.set_ylabel(chrom, fontsize=10, fontweight='bold')
        ax.set_yticks([])
        
        # Remove x-axis ticks and labels (no grid, no axis clutter)
        ax.set_xticks([])
        ax.set_xticklabels([])
        
        # Remove spines (like haplopainting does)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
    # Add x-axis only at the bottom with round number marks
    # Generate round number ticks
    range_mbp = (global_max - global_min) / 1e6
    if range_mbp <= 100:
        tick_interval = 10
    elif range_mbp <= 300:
        tick_interval = 20
    elif range_mbp <= 700:
        tick_interval = 50
    else:
        tick_interval = 100
    
    # Create tick positions at round numbers
    first_tick = int((global_min / 1e6 + tick_interval) // tick_interval) * tick_interval
    ticks = []
    labels = []
    tick_pos = first_tick
    while tick_pos * 1e6 <= global_max:
        ticks.append(tick_pos * 1e6)
        labels.append(str(tick_pos))
        tick_pos += tick_interval
    
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels(labels)
    axes[-1].spines['bottom'].set_visible(True)
    axes[-1].set_xlabel("Position (Mbp)", fontsize=12)
    
    # Add overall title
    fig.suptitle(f"Imputed - {base_name}", fontsize=14, fontweight='bold', y=0.98)
    
    # Create legend with unique genomes
    unique_genomes = sorted(dataframe['genome'].unique())
    legend_patches = [mpatches.Patch(facecolor=colors_dict[genome], 
                                     edgecolor='none', label=genome)
                     for genome in unique_genomes]
    
    fig.legend(handles=legend_patches, loc='upper right', fontsize=8, ncol=1, framealpha=0.95)
    
    # Save with minimal spacing - very tight spacing between chromosomes
    output_file = os.path.join(output_folder, f"{base_name}.png")
    plt.subplots_adjust(left=0.08, right=0.85, top=0.96, bottom=0.01, hspace=0.01)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    print(f"Plot saved: {output_file}")


def create_summary_statistics(dataframe, verbose=False):
    """
    Calculate and print summary statistics if verbose mode is enabled.
    
    Args:
        dataframe (pd.DataFrame): DataFrame with haplotype ranges
        verbose (bool): If True, print statistics to terminal
    """
    if not verbose:
        return
    
    # Calculate total BP across all genomes for percentage calculation
    total_bp_all = (dataframe['ref_end'] - dataframe['ref_start']).sum()
    
    # Calculate statistics per genome
    stats = []
    for genome in sorted(dataframe['genome'].unique()):
        genome_data = dataframe[dataframe['genome'] == genome]
        total_bp = (genome_data['ref_end'] - genome_data['ref_start']).sum()
        num_ranges = len(genome_data)
        pct_coverage = (total_bp / total_bp_all * 100) if total_bp_all > 0 else 0
        
        stats.append({
            'Genome': genome,
            'Total_BP_Covered': total_bp,
            'Num_Ranges': num_ranges,
            'Percent_Of_Total': f"{pct_coverage:.2f}%"
        })
    
    stats_df = pd.DataFrame(stats)
    print("\nSummary Statistics:")
    print(stats_df.to_string(index=False))




def main(args=None):
    """
    Main function to handle argument parsing and execution flow.
    
    Args:
        args: List of command line arguments (for testing/module integration)
    """
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="plot-imputed-hvcf",
        description="Plot imputed hVCF files showing genome-colored haplotype ranges across all chromosomes."
    )

    parser.add_argument("pangenome_hvcf_folder", help="Path to folder containing pangenome hVCF files.")
    parser.add_argument("imputed_hvcf", help="Path to the imputed hVCF file.")
    parser.add_argument("reference_hvcf", help="Path to the reference hVCF file.")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output directory for plots (default: pangenome_folder/plots)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print detailed processing messages.")

    parsed_args = parser.parse_args(args)

    # Set plot folder
    plot_folder = parsed_args.output if parsed_args.output else f"{parsed_args.pangenome_hvcf_folder}/plots/"
    
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

    # Get list of genomes from pangenome folder
    hvcf_files_path = [os.path.join(parsed_args.pangenome_hvcf_folder, file) 
                       for file in os.listdir(parsed_args.pangenome_hvcf_folder) 
                       if file.endswith('.h.vcf.gz')]
    genomes = [os.path.basename(file).split('.')[0] for file in hvcf_files_path]
    reference_name = os.path.basename(parsed_args.reference_hvcf).split('.')[0]

    if parsed_args.verbose:
        print(f"\nImputed hVCF file: {parsed_args.imputed_hvcf}")
        print(f"Pangenome genomes ({len(genomes)}): {genomes}")
        print(f"Reference genome: {reference_name}")
        print(f"Output folder: {plot_folder}\n")

    # Assign colors
    colors_dict = assign_colors(genomes, reference_name)

    # Parse hVCF to dataframe
    if parsed_args.verbose:
        print("Parsing hVCF file...")
    dataframe = parse_hvcf_to_dataframe(parsed_args.imputed_hvcf, verbose=parsed_args.verbose)
    
    if dataframe.empty:
        print("ERROR: No haplotype ranges found in hVCF file!")
        sys.exit(1)
    
    if parsed_args.verbose:
        print(f"Found {len(dataframe)} haplotype ranges across {len(dataframe['chr'].unique())} chromosomes")

    # Create plot
    if parsed_args.verbose:
        print("Creating haplotype painting plot...")
    plot_haplotype_painting_with_data(dataframe, colors_dict, plot_folder, parsed_args.imputed_hvcf, verbose=parsed_args.verbose)

    # Create statistics (verbose only)
    create_summary_statistics(dataframe, verbose=parsed_args.verbose)

    print(f"\nPlotting complete!")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nScript interrupted by user. Exiting...")
        sys.exit(0)
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
