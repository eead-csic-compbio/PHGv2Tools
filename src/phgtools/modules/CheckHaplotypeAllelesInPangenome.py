#!/usr/bin/env python
"""
A module to query a pangenome hapIDranges.tsv file for lines that overlap
a specified coordinate and prints the full overlapping lines, including the header.
"""

import sys
import argparse
import gzip

def DefineRange(hapid_ranges_file, start, end, chromosome, verbose=False):
    """
    Identifies and prints the lines (genomic ranges) in the hapIDranges.tsv
    that overlap the user-specified coordinates.
    
    Args:
        hapid_ranges_file (str): Path to the merged hapIDranges.tsv file.
        start (int): Start coordinate.
        end (int): End coordinate.
        chromosome (str): Chromosome number (e.g., "3" for "chr3").
        verbose (bool): If True, print debug statements for checks.
    """
    
    # Prepend 'chr' for standard VCF format (TSV typically uses chrX, e.g., chr1H)
    chromosome_label = chromosome
    
    if verbose:
        print(f"Querying for ranges on {chromosome_label}: {start} to {end}\n")
        print("--- VERBOSE MODE: Checking TSV records ---")

    overlapping_lines = []
    header_line = None
    line_number = 0

    # Auto-detect gzipped files
    if hapid_ranges_file.endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt')  # 'rt' = read text mode
    else:
        open_func = lambda f: open(f, 'r')

    with open_func(hapid_ranges_file) as f:
        for line in f:
            line_number += 1
            line = line.strip()
            
            if line.startswith("#"):
                # Capture the header (starts with #CHROM) or any other initial comment lines
                header_line = line
                if line.startswith("#CHROM"):
                    print(line) # Print the header immediately
                continue
            
            fields = line.split("\t")
            
            # TSV fields: CHROM=fields[0], START=fields[1], END=fields[2]
            
            current_chrom = fields[0]
            
            # 1. Check if the line is on the correct chromosome                
            if not current_chrom == chromosome_label:
                continue
            
            try:
                ref_start = int(fields[1])
                ref_end = int(fields[2])
            except ValueError:
                if verbose:
                    print(f"[{line_number}] WARNING: Could not parse START/END coordinates. Skipping.")
                continue
                            
            # Check for overlap: max(ref_start, start) < min(ref_end, end)
            is_overlap = (ref_start <= end and ref_end >= start) 
            
            if is_overlap:                
                # Print the overlapping line to the terminal
                print(line) 


def main(args=None):
    """
    Main function to handle script execution and user input/arguments.
    
    Args:
        args: List of command line arguments (for testing/module integration)
    """
    
    # Use provided args or sys.argv
    if args is None:
        args = sys.argv[1:]
    
    # Argument Parsing
    parser = argparse.ArgumentParser(
        prog="check-haplotype-alleles",
        description="Query a pangenome hapIDranges.tsv to find overlapping ranges and print the full lines." \
        "\n\nIf you do not provide coordinates via command line, you will be prompted to enter them.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments: Only need the TSV file now
    parser.add_argument("hapid_ranges_file", help="Path to the merged hapIDranges.tsv file.")
    
    # Optional arguments
    parser.add_argument("-s", "--start", type=int, help="Start coordinate (e.g., 10000).")
    parser.add_argument("-e", "--end", type=int, help="End coordinate (e.g., 10500).")
    parser.add_argument("-c", "--chromosome", help="Chromosome name, exactly how it appears in the hapIDranges.tsv file (e.g., 'chr1H').")
    
    # Verbose Flag
    parser.add_argument("-v", "--verbose", action="store_true", help="Print debug statements showing all checks.")
    
    parsed_args = parser.parse_args(args)

    # Coordinate Input and Argument Logging
    
    # Use the first 3 positional arguments for the coordinate input
    start, end, chromosome, verbose = parsed_args.start, parsed_args.end, parsed_args.chromosome, parsed_args.verbose
    
    if start is None or end is None or chromosome is None:
        print("Coordinates not fully provided via command line. Please enter them now:")
        try:
            if start is None:
                start = int(input("Start: "))
            if end is None:
                end = int(input("End: "))
            if chromosome is None:
                # Assuming the user enters "1H" for chr1H, etc.
                chromosome = input("Chromosome (enter only the number/label, e.g., '1H'): ")
        except ValueError:
            print("Error: Coordinates must be valid integers.")
            sys.exit(1)
        except KeyboardInterrupt:
            print("\nScript interrupted by user. Exiting...")
            sys.exit(0)

    # Verbose: Print arguments used
    if verbose:
        print("\n--- SCRIPT ARGUMENTS ---")
        print(f"HapID Ranges File:  {parsed_args.hapid_ranges_file}")
        print(f"Query Chromosome:   {chromosome}")
        print(f"Query Start:        {start}")
        print(f"Query End:          {end}")
        print("------------------------\n")

    # Core Logic Execution
    
    # The function now handles the printing directly
    DefineRange(parsed_args.hapid_ranges_file, start, end, chromosome, verbose)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nScript interrupted by user. Exiting...")
        sys.exit(0)
    except FileNotFoundError:
        print("\nError: One or more input files were not found. Please check paths.")
        sys.exit(1)
