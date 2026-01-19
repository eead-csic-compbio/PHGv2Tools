#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
vcfDistance.py
Generate distance matrix comparing all varieties in a g.VCF file.
Outputs distance matrix (1 = identity) and count matrix (number of common SNPs).
Includes optional heatmap visualization with hierarchical clustering.
"""

import re
import subprocess
import argparse
import os
import gzip
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import numpy as np
import warnings
from multiprocessing import Pool, cpu_count


def open_vcf_file(vcf_file):
    """Opens VCF file handling both bgzip compressed and plain text files."""
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')


def get_variety_columns(VCFfile, verbose=False):
    """Reads VCF header to get genotype column indices for all variety."""
    VCFcolumns = {}
    
    # Open the file (handles both bgzip and plain text)
    try:
        with open_vcf_file(VCFfile) as file:
            for line in file:
                if line.startswith('#CHROM\t'):
                    fields = line.split()
                    # Create dict {variety_name: column_index} starting after the 9 fixed columns (i + 10)
                    # The VCF format is 1-based, and the samples start at column 10 (index 9)
                    VCFcolumns = {field: i + 10 for i, field in enumerate(fields[9:])}
                    break
    except Exception as e:
        sys.exit(f"Error reading VCF header: {e}")

    if not VCFcolumns:
        sys.exit("Error: Could not parse VCF header or find genotype columns.")
        
    if verbose:
        print("\n--- Variety Selection ---")
        print("Processing ALL varieties present in the VCF file (default):")
        for key in VCFcolumns.keys():
            print(f"- {key}")
        print("-------------------------\n")
        
    return VCFcolumns


def calculate_pair_distance(args):
    """
    Worker function for parallel processing of variety pairs.
    Calculates distance and count for a single pair comparison.
    
    Args:
        args: tuple of (vcf_file, col_b1, col_b2, b1_name, b2_name)
    
    Returns:
        tuple: (b1_name, b2_name, eq, neq)
    """
    vcf_file, col_b1, col_b2, b1_name, b2_name = args
    
    eq = 0
    neq = 0
    
    # Use appropriate decompression based on file type
    if vcf_file.endswith('.gz'):
        cmd = f"zcat {vcf_file} | cut -f {col_b1},{col_b2}"
    else:
        cmd = f"cat {vcf_file} | cut -f {col_b1},{col_b2}"
    
    try:
        # Execute command and read output line by line
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        for line in process.stdout:
            if line.startswith('#'): continue
                
            fields = line.split()

            if len(fields) < 2: continue
                
            gt1 = fields[0].strip()
            gt2 = fields[1].strip()

            # Skip missing data ('.')
            if gt1 == '.' or gt2 == '.': continue
            
            try:
                gt1 = int(gt1)
                gt2 = int(gt2)
            except ValueError:
                
                # Skip if any of the GT values are heterozygous or contain missing data characters (., /, |)
                if '/' in gt1 or '|' in gt1 or '.' in gt1 or '/' in gt2 or '|' in gt2 or '.' in gt2:
                    continue
                
                gt1 = gt1.split(':')[0]
                gt2 = gt2.split(':')[0]
                if len(gt1) != 1 or len(gt2) != 1: continue # Skip if not a single number (e.g. 0/0, 1|1)
                
                # Final check that we can compare the alleles (e.g., '0' vs '1')
                if gt1 == gt2:
                    eq += 1
                else:
                    neq += 1
                continue
                
            # If the GT fields were simple integers
            if gt1 == gt2:
                eq += 1
            else:
                neq += 1
                
        # Check for external process errors
        process.stdout.close()
        process.wait()
        if process.returncode != 0:
            stderr_output = process.stderr.read()
            print(f"Error executing subprocess command: {cmd}", file=sys.stderr)
            print(f"Subprocess STDERR: {stderr_output}", file=sys.stderr)
            return (b1_name, b2_name, None, None)

    except subprocess.CalledProcessError as e:
        print(f"Error during VCF processing: {e}", file=sys.stderr)
        return (b1_name, b2_name, None, None)
    
    return (b1_name, b2_name, eq, neq)


def calculate_distance_matrices(VCFfile, VCFcolumns, distfile, countfile, verbose=False, num_threads=1):
    """
    Calculates the distance and count matrices by streaming data from the VCF.
    Uses parallel processing to compare multiple variety pairs simultaneously.
    
    Args:
        VCFfile: Path to input VCF file
        VCFcolumns: Dictionary of {variety_name: column_index}
        distfile: Output distance matrix file path
        countfile: Output count matrix file path
        verbose: Print processing messages
        num_threads: Number of threads to use (default: 1 thread)
    """
    if num_threads is None or num_threads <= 0:
        num_threads = 1
    
    variety = sorted(VCFcolumns.keys())
    
    # Open output files
    try:
        distmat = open(distfile, 'w')
        countmat = open(countfile, 'w')
    except IOError as e:
        sys.exit(f"Error opening output files: {e}")

    # Write the header line for both matrices
    header = "\t" + "\t".join(variety) + "\n"
    distmat.write(header)
    countmat.write(header)

    if verbose:
        print(f"\n--- Distance Matrix Calculation (using {num_threads} threads) ---")
    
    # Prepare all comparison pairs to process
    # We only need to compute the upper triangle (symmetric matrix)
    comparison_pairs = []
    for i, b1 in enumerate(variety):
        for j, b2 in enumerate(variety):
            if i >= j:  # Skip diagonal and lower triangle (will fill from cached results)
                continue
            col_b1 = VCFcolumns[b1]
            col_b2 = VCFcolumns[b2]
            comparison_pairs.append((VCFfile, col_b1, col_b2, b1, b2))
    
    # Run parallel processing
    results_dict = {}
    with Pool(num_threads) as pool:
        if verbose:
            print(f"Processing {len(comparison_pairs)} variety pairs in parallel...")
        
        results = pool.imap_unordered(calculate_pair_distance, comparison_pairs)
        for b1, b2, eq, neq in results:
            if eq is not None:
                # Store result for both (b1,b2) and (b2,b1) since matrix is symmetric
                results_dict[(b1, b2)] = (eq, neq)
                results_dict[(b2, b1)] = (eq, neq)
                if verbose:
                    print(f"  ✓ {b1} vs {b2}: {eq} matches, {neq} mismatches")
    
    # Write output matrices with results in order
    for i, b1 in enumerate(variety):
        row_output_dist = [b1]
        row_output_count = [b1]
        
        for j, b2 in enumerate(variety):
            if b1 == b2:
                row_output_dist.append("0.0000")
                row_output_count.append("NA")
            else:
                # Get cached result
                eq, neq = results_dict.get((b1, b2), (0, 0))
                
                # Output the distance and count
                if eq + neq == 0:
                    distance = "NA"
                    count = "0"
                else:
                    distance = 1 - (eq / (eq + neq))
                    count = eq + neq
                    # Format distance to 4 decimal places for output
                    distance = f"{distance:.4f}"
                
                row_output_dist.append(distance)
                row_output_count.append(str(count))
        
        # Write the completed row to the files
        distmat.write("\t".join(row_output_dist) + "\n")
        countmat.write("\t".join(row_output_count) + "\n")
        distmat.flush()
        countmat.flush()
        
    distmat.close()
    countmat.close()
    
    if verbose:
        print(f"\nAll done! Check your files: {distfile} and {countfile}")


def plot_heatmap(distfile, plotfile, verbose=False):
    """
    Generates a clustered heatmap with dendrogram on the left side using seaborn.clustermap.
    This ensures perfect alignment between dendrogram and heatmap rows.
    """
    
    if verbose:
        print(f"\n--- Creating Clustered Heatmap with Dendrogram ---")
        print(f"Saving to: {plotfile}")

    try:
        # Read the generated distance matrix
        df_dist = pd.read_csv(distfile, sep="\t", index_col=0)
        
        # Fix for header misalignment if it happens
        if df_dist.columns[0].startswith('Unnamed'):
             df_dist = pd.read_csv(distfile, sep="\t", index_col=0, header=0)
             df_dist.columns = df_dist.index
        
        # Convert to float and handle 'NA'
        df_dist = df_dist.replace('NA', np.nan)
        
        # Ensure all columns are numeric before plotting
        for col in df_dist.columns:
            df_dist[col] = pd.to_numeric(df_dist[col], errors='coerce')
        
        # Fill diagonal (which should be 0) just in case
        np.fill_diagonal(df_dist.values, 0)
        
        # Calculate color scale range from non-diagonal values only (exclude 0s on diagonal)
        mask_diag = np.eye(len(df_dist), dtype=bool)
        off_diag_values = df_dist.values[~mask_diag]
        vmin = np.nanmin(off_diag_values)
        vmax = np.nanmax(off_diag_values)
        
    except Exception as e:
        sys.exit(f"Error reading distance matrix for plotting: {e}")

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        
        # Calculate number of samples for dynamic sizing
        n_samples = len(df_dist)
        cell_size = 0.6  
        fig_width = max(10, n_samples * cell_size + 2)  # Extra space for dendrogram
        fig_height = max(8, n_samples * cell_size)
        
        # Create custom annotation array - empty string for diagonal
        annot_labels = df_dist.applymap(lambda x: f'{x:.2f}')
        for i in range(len(df_dist)):
            annot_labels.iloc[i, i] = ''  # No text on diagonal
        
        # Create a continuous colormap (blue = similar, red = distant)
        from matplotlib.colors import LinearSegmentedColormap
        colors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', 
                  '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
        cmap_smooth = LinearSegmentedColormap.from_list('smooth_diverging', colors, N=256)
        
        # Create clustermap with BOTH row and column clustering
        # This ensures the diagonal stays as a diagonal after reordering
        g = sns.clustermap(
            df_dist,
            cmap=cmap_smooth,          # Smooth continuous diverging colormap
            annot=annot_labels,        # Custom annotations (no 0.00 on diagonal)
            fmt='',                    # No formatting needed, already formatted
            vmin=vmin,                 # Color scale min (excluding diagonal)
            vmax=vmax,                 # Color scale max (excluding diagonal)
            figsize=(fig_width, fig_height),
            dendrogram_ratio=(0.15, 0.08),  # Left dendrogram larger, small top
            row_cluster=True,          # Cluster rows
            col_cluster=True,          # Cluster columns (keeps diagonal aligned)
            cbar_pos=None,             # No colorbar
            linewidths=0.5,            # Cell borders
            linecolor='#f0f0f0',       # Light gray cell borders for subtle separation
            tree_kws={
                'linewidths': 2.5,     # Thicker dendrogram lines
                'colors': None,        # Will be set by color_threshold
            },
            annot_kws={
                'size': 9,             # Annotation font size
                'weight': 'semibold',
                'color': '#333333'     # Dark gray text for better readability
            }
        )
        
        # Paint diagonal cells black (self-comparisons = 0.00)
        # The diagonal in the reordered matrix is still at position (i, i)
        for i in range(n_samples):
            g.ax_heatmap.add_patch(plt.Rectangle(
                (i, i), 1, 1,  # x, y, width, height - diagonal is (i, i)
                fill=True, facecolor='black', edgecolor='white', linewidth=0.3
            ))
        
        # Hide the top dendrogram (we only want left side visible)
        g.ax_col_dendrogram.set_visible(False)
        
        # Style the dendrogram axis
        g.ax_row_dendrogram.set_xticks([])
        g.ax_row_dendrogram.set_yticks([])
        for spine in g.ax_row_dendrogram.spines.values():
            spine.set_visible(False)
        
        # Remove axis labels (variety names appear on axes)
        g.ax_heatmap.set_xlabel("")
        g.ax_heatmap.set_ylabel("")
        
        # Improve tick label appearance
        g.ax_heatmap.tick_params(axis='both', which='major', labelsize=10)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontweight='medium')
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontweight='medium')
        
        # Adjust layout for better spacing
        g.fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.15)

    # Save as PNG at 300 DPI
    g.savefig(plotfile, dpi=300, bbox_inches='tight', format='png')
    plt.close(g.fig)

    if verbose:
        print(f"Visualization saved to: {plotfile}")


def main(args=None):
    """
    Main function to handle argument parsing and execution flow.
    
    Args:
        args: List of command line arguments (for testing/module integration)
    """
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="vcf-distance",
        description="Generate a distance matrix comparing all varieties in a g.VCF file. "
                    "The g.vcf file should be an output of a pangenome pipeline (e.g., PHG) "
                    "and can be merged and bgzip compressed (.vcf.gz) or plain text (.vcf). "
                    "Outputs distance matrix (1 = identity) and count matrix (number of common SNPs).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\nHeatmap and Dendrogram Visualization:\n"
               "  Use -p/--heatmap-plot to generate a publication-quality PNG (300 DPI) with:\n"
               "    - Left side: Hierarchical Clustering Dendrogram showing phylogenetic relationships\n"
               "    - Right side: Distance Heatmap (reordered) with variety names aligned to dendrogram\n"
               "    - X and Y axes both use dendrogram-ordered variety names for perfect alignment\n"
               "\nParallel Processing:\n"
               "  Use -t/--threads to speed up processing with multiple CPU cores:\n"
               "    - Default: 1 thread (sequential processing)\n"
               "    - Specify number: -t 4 to use 4 cores\n"
               "    - Use all cores: -t -1\n"
               "\nExample usage:\n"
               "  phgtools vcf-distance input.vcf.gz -p              # Generate with heatmap (1 thread)\n"
               "  phgtools vcf-distance input.vcf.gz -t 8 -v         # Use 8 threads with verbose output\n"
               "  phgtools vcf-distance input.vcf.gz -t 4 -p         # 4 threads + heatmap visualization"
    )
    
    parser.add_argument("vcf_file", help="Input VCF file (.vcf.gz compressed or plain .vcf text)")
    parser.add_argument("-o", "--out-matrix", type=str, default="distance_matrix.tsv",
                        help="Output distance matrix file path (default: distance_matrix.tsv). "
                             "Count matrix will be named with '_count' suffix before .tsv")
    parser.add_argument("-p", "--heatmap-plot", nargs='?', const='distance_heatmap.png', default=None, type=str,
                        help="Generate heatmap with dendrogram visualization (PNG, 300 DPI). "
                             "If no path provided, uses 'distance_heatmap.png'. "
                             "Can also provide existing distance matrix .tsv file to re-plot.")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for parallel processing (default: 1 thread). "
                             "Use -1 to use all available cores, or specify a number like -t 4.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print detailed processing messages")

    parsed_args = parser.parse_args(args)

    # Validate input file
    vcf_file = parsed_args.vcf_file
    if not os.path.exists(vcf_file):
        sys.exit(f"Error: VCF file not found: {vcf_file}\n"
                 f"Supported formats: .vcf.gz (bgzip compressed) or .vcf (plain text)")

    distfile = parsed_args.out_matrix
    # Generate count matrix name by inserting '_count' before .tsv extension
    if distfile.endswith('.tsv'):
        countfile = distfile.replace('.tsv', '_count.tsv')
    else:
        countfile = distfile + '_count.tsv'

    # Handle thread count
    num_threads = parsed_args.threads
    if num_threads == -1:
        num_threads = cpu_count()
    elif num_threads is None or num_threads <= 0:
        num_threads = 1
    
    if parsed_args.verbose:
        print(f"\nInput VCF file: {vcf_file}")
        print(f"Output distance matrix: {distfile}")
        print(f"Output count matrix: {countfile}")
        print(f"Using {num_threads} thread(s) for processing\n")

    # Check if user provided an existing matrix file to plot
    if parsed_args.heatmap_plot is not None and parsed_args.heatmap_plot != 'distance_heatmap.png' and parsed_args.heatmap_plot.endswith('.tsv'):
        # User provided an existing matrix file to plot
        if not os.path.exists(parsed_args.heatmap_plot):
            sys.exit(f"Error: The provided distance matrix file does not exist: {parsed_args.heatmap_plot}")
        
        if parsed_args.verbose:
            print(f"Using existing distance matrix file: {parsed_args.heatmap_plot}")
        
        # Use provided matrix for plotting
        plot_heatmap(parsed_args.heatmap_plot, 'distance_heatmap.png', verbose=parsed_args.verbose)
        
    else:
        # Normal workflow: generate matrices and optionally plot
        
        # Get variety columns from VCF header
        VCFcolumns = get_variety_columns(vcf_file, verbose=parsed_args.verbose)
        
        # Calculate distance matrices with parallel processing
        calculate_distance_matrices(vcf_file, VCFcolumns, distfile, countfile, verbose=parsed_args.verbose, num_threads=num_threads)
        
        # Plot if requested
        if parsed_args.heatmap_plot is not None:
            plot_heatmap(distfile, parsed_args.heatmap_plot, verbose=parsed_args.verbose)

    if parsed_args.verbose:
        print("\nDistance matrix analysis complete!")


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
