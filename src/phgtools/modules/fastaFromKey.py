import os
import argparse
import gzip
import sys
import subprocess
import shutil

def extract_coord_from_key(vcf_file, key, vcf_folder, verbose=False):
    """
    Parses the input file (VCF file or VCF folder) to find the coordinates
    associated with the given haplotype key.
    """
    sample_name = None
    chr_ = None
    start = None
    end = None
    found = False

    # Option A: Search in VCF file
    if vcf_file is not None:
        if verbose:
            print(f"# Searching for key {key} in VCF file: {vcf_file}")
        try:
            # Handle both gzip and plain text for vcf
            open_func = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with open_func(vcf_file, mode) as f:
                for line in f:
                    if line.startswith("##ALT"):
                        if key in line:
                            # Parse ##ALT line
                            # Example: ##ALT=<ID=hash,Description="...",SampleName=Name,Regions=chr:s-e,...>
                            parts = line.split(",")
                            for part in parts:
                                if part.strip().startswith("SampleName="):
                                    sample_name = part.split("=")[1]
                                elif part.strip().startswith("Regions="):
                                    region = part.split("=")[1].strip('"')
                                    # Handle standard region format chr:start-end
                                    if ":" in region and "-" in region:
                                        chr_ = region.split(":")[0]
                                        coords = region.split(":")[1]
                                        start = coords.split("-")[0]
                                        end = coords.split("-")[1]
                            found = True
                            break
        except Exception as e:
            print(f"Error reading VCF file: {e}", file=sys.stderr)
            sys.exit(1)

    # Option B: Search in VCF folder (file by file) if VCF file not provided
    elif vcf_folder is not None and not found:
        if verbose:
            print(f"# Searching for key {key} in VCF files within: {vcf_folder}")
        try:
            if not os.path.exists(vcf_folder):
                print(f"Error: VCF folder not found: {vcf_folder}", file=sys.stderr)
                sys.exit(1)
            
            vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith(('.h.vcf', '.h.vcf.gz'))]
            
            if not vcf_files:
                print(f"Error: No VCF files found in {vcf_folder}", file=sys.stderr)
                sys.exit(1)
            
            for vcf_filename in vcf_files:
                if found:
                    break
                
                vcf_path = os.path.join(vcf_folder, vcf_filename)
                open_func = gzip.open if vcf_filename.endswith('.gz') else open
                mode = 'rt' if vcf_filename.endswith('.gz') else 'r'
                
                with open_func(vcf_path, mode) as f:
                    for line in f:
                        if line.startswith("##ALT"):
                            if key in line:
                                # Parse ##ALT line from VCF style
                                parts = line.split(",")
                                for part in parts:
                                    if part.strip().startswith("SampleName="):
                                        sample_name = part.split("=")[1]
                                    elif part.strip().startswith("Regions="):
                                        region = part.split("=")[1].strip('"')
                                        if ":" in region and "-" in region:
                                            chr_ = region.split(":")[0]
                                            coords = region.split(":")[1]
                                            start = coords.split("-")[0]
                                            end = coords.split("-")[1]
                                found = True
                                break
                        # Also check in VCF body (ALT column contains the hash)
                        elif not line.startswith("#"):
                            cols = line.split("\t")
                            if len(cols) > 7:
                                alt_val = cols[4]
                                if key in alt_val:
                                    sample_name = vcf_filename.replace('.h.vcf.gz', '').replace('.h.vcf', '')
                                    chr_ = cols[0]
                                    start = cols[1]
                                    # INFO field (col 7) usually contains END=xxx
                                    info_field = cols[7]
                                    end_match = [x for x in info_field.split(";") if x.startswith("END=")]
                                    if end_match:
                                        end = end_match[0].split("=")[1]
                                    found = True
                                    break
            
            if found and verbose:
                print(f"# Found in file: {vcf_filename}")
        except Exception as e:
            print(f"Error searching VCF folder: {e}", file=sys.stderr)
            sys.exit(1)

    if not found:
        print(f"Error: Key {key} not found in the provided files.", file=sys.stderr)
        sys.exit(1)

    if sample_name and chr_ and start and end:
        if verbose:
            print(f"# Found -> Sample: {sample_name} | Location: {chr_}:{start}-{end}")
        return sample_name, chr_, start, end
    else:
        print("Error: Could not parse coordinates from the found entry.", file=sys.stderr)
        sys.exit(1)


def extract_fasta_from_coord(fastas_folder, sample_name, chr_, start, end, output_folder, key, verbose=False):
    """
    Uses samtools faidx to extract the sequence.
    """
    # Check for samtools
    if shutil.which("samtools") is None:
        print("Error: 'samtools' is not installed or not in your PATH.", file=sys.stderr)
        sys.exit(1)

    fasta_file = os.path.join(fastas_folder, f"{sample_name}.fa")
    
    if not os.path.exists(fasta_file):
        # Try checking for .fasta extension as fallback
        fasta_file_alt = os.path.join(fastas_folder, f"{sample_name}.fasta")
        if os.path.exists(fasta_file_alt):
            fasta_file = fasta_file_alt
        else:
            print(f"Error: Fasta file not found: {fasta_file}", file=sys.stderr)
            sys.exit(1)


    # Convert start and end to integers for comparison
    try:
        start_int = int(start)
        end_int = int(end)
    except Exception:
        print(f"Error: Could not convert start/end to integer: start={start}, end={end}", file=sys.stderr)
        sys.exit(1)

    # If start > end, print warning and swap
    if start_int > end_int:
        print(f"Warning: start ({start_int}) > end ({end_int}), swapping positions.", file=sys.stderr)
        start_int, end_int = end_int, start_int
    region_string = f"{chr_}:{start_int}-{end_int}"
    if verbose:
        print(f"\n# Extracting sequence for {sample_name} at {region_string}...")

    cmd = ["samtools", "faidx", fasta_file, region_string]

    if output_folder:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        
        output_filename = f"{key}.fa"
        output_path = os.path.join(output_folder, output_filename)
        
        try:
            # Get the sequence from samtools
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse the output to modify the header
            lines = result.stdout.strip().split('\n')
            with open(output_path, "w") as outfile:
                # Write custom header
                outfile.write(f">{key}_{sample_name}@{chr_}:{start}-{end}\n")
                # Write sequence (skip the original header)
                for line in lines[1:]:
                    outfile.write(line + "\n")
            
            if verbose:
                print(f"# Sequence saved to: {output_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools: {e}", file=sys.stderr)
            sys.exit(1)
        except IOError as e:
            print(f"Error writing output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Print to stdout with modified header
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            lines = result.stdout.strip().split('\n')
            # Print custom header
            print(f">{key}_{sample_name}@{chr_}:{start}-{end}")
            # Print sequence (skip the original header)
            for line in lines[1:]:
                print(line)
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools: {e}", file=sys.stderr)
            sys.exit(1)


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Extract FASTA sequence from a PHGv2 key.",
        epilog="Requires 'samtools' to be installed and available in PATH."
    )

    # Required Arguments
    parser.add_argument("--key", required=True, help="The haplotype key/hash to extract.")
    parser.add_argument("--fastas-folder", required=True, help="Folder containing the genome assemblies (.fa files).")

    # Mutually Exclusive: Need either VCF file or VCF folder
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vcf-file", help="Path to a VCF file to search in.")
    group.add_argument("--vcf-folder", help="Folder containing individual h.vcf files to search through.")

    # Optional
    parser.add_argument("--output-folder", help="Folder to save the output FASTA. If not provided, prints to terminal.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output.")

    parsed_args = parser.parse_args(args)

    # Execute Logic
    sample_name, chr_, start, end = extract_coord_from_key(
        parsed_args.vcf_file,
        parsed_args.key,
        parsed_args.vcf_folder,
        verbose=parsed_args.verbose
    )

    extract_fasta_from_coord(
        parsed_args.fastas_folder,
        sample_name,
        chr_,
        start,
        end,
        parsed_args.output_folder,
        parsed_args.key,
        verbose=parsed_args.verbose
    )

if __name__ == "__main__":
    main()
