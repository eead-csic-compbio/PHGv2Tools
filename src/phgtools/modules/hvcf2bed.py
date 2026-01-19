import sys
import gzip
import os
import re
import argparse
import glob

def parse_vcf_lines(lines, genome_name, verbose=False):
    viewed_checksums = set() 
    empty_ranges = 0
    skipped_checksums_bug = 0
    output_lines = []
    prev_checksum = None

    # --- REGEX PATTERNS ---
    
    # 1. Multi-region format (quoted Regions)
    pat_multi = re.compile(
        r'SampleName=([^,>]+).*?Regions="([^"]+)".*?Checksum=([^,]+).*?RefChecksum=([^,>]+).*?RefRange=([^:]+):(\d+)-(\d+)'
    )
    
    # 2. Standard format (RefChecksum BEFORE RefRange)
    pat_std = re.compile(
        r'SampleName=([^,]+).*?Regions=([^:]+):(\d+)-(\d+)(?!.*Regions=).*?Checksum=([^,]+).*?RefChecksum=([^,]+).*?RefRange=([^:]+):(\d+)-(\d+)'
    )

    # 3. SWAPPED format (RefRange BEFORE RefChecksum) - This fixes your Ler0_1 issue
    pat_swapped = re.compile(
        r'SampleName=([^,]+).*?Regions=([^:]+):(\d+)-(\d+)(?!.*Regions=).*?Checksum=([^,]+).*?RefRange=([^:]+):(\d+)-(\d+).*?RefChecksum=([^,>]+)'
    )

    def process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum, is_multi=False):
        nonlocal prev_checksum, empty_ranges, skipped_checksums_bug, verbose
        
        # Validate RefRange logic
        try:
            ref_start, ref_end = int(ref_start), int(ref_end)
        except ValueError:
            return

        # For multi-region entries, we rely on RefRange being valid
        if is_multi:
            ref_length = abs(ref_end - ref_start)
            if ref_length <= 0:
                empty_ranges += 1
                return
            
            output_lines.append(
                f"{ref_chr}\t{ref_start}\t{ref_end}\t+\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
            )
            return
        
        # For standard/swapped formats
        start, end = int(start), int(end)
        length = end - start
        
        strand = "+"
        if start > end:
            strand = "-"
            start, end = end, start
            length = end - start
        
        if length < 0:
            empty_ranges += 1
            return 

        # Check for duplicates
        if checksum == prev_checksum:
            skipped_checksums_bug += 1
            return 
        elif checksum in viewed_checksums:
            pass # Keep going if seen before but not sequential

        output_lines.append(
            f"{chr_}\t{start}\t{end}\t{strand}\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
        )
        
        prev_checksum = checksum
        viewed_checksums.add(checksum)

    # --- MAIN PARSING LOOP ---
    for lineno, line in enumerate(lines, start=1):
        line = line.strip()
        
        if not line.startswith("##ALT"):
            continue

        # Try Pattern 1: Multi-region
        m2 = pat_multi.search(line)
        if m2:
            genome, regions, checksum, ref_checksum, ref_chr, ref_start, ref_end = m2.groups()
            process_entry(None, None, None, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum, is_multi=True)
            prev_checksum = checksum
            continue

        # Try Pattern 2: Standard (RefChecksum first)
        m1 = pat_std.search(line)
        if m1:
            genome, chr_, start, end, checksum, ref_checksum, ref_chr, ref_start, ref_end = m1.groups()
            process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum)
            continue

        # Try Pattern 3: Swapped (RefRange first)
        m3 = pat_swapped.search(line)
        if m3:
            genome, chr_, start, end, checksum, ref_chr, ref_start, ref_end, ref_checksum = m3.groups()
            process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum)
            continue

        if verbose:
            print(f"# WARNING (line {lineno}): No regex match for ALT line.")

    if verbose:
        print(f"# Processed {len(output_lines)} haplotype blocks")

    return output_lines

def main(args=None):
    parser = argparse.ArgumentParser(description='Convert HVCF files to BED format')
    parser.add_argument('vcf_folder', help='Path to folder containing h.vcf files')
    parser.add_argument('genome_name', nargs='?', default=None, help='Genome name (optional)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    vcf_folder = args.vcf_folder
    genome_name = args.genome_name
    verbose = getattr(args, 'verbose', False)

    if not os.path.isdir(vcf_folder):
        print(f"ERROR: Folder not found: {vcf_folder}")
        sys.exit(1)

    if genome_name is None:
        # Find all files
        files_gz = glob.glob(os.path.join(vcf_folder, "*.h.vcf.gz"))
        files_plain = glob.glob(os.path.join(vcf_folder, "*.h.vcf"))
        
        # Unique base names
        targets = set()
        for f in files_gz + files_plain:
            base = os.path.basename(f).replace('.h.vcf.gz', '').replace('.h.vcf', '')
            targets.add(base)
            
        for target in targets:
            process_single_file(vcf_folder, target, verbose)
    else:
        process_single_file(vcf_folder, genome_name, verbose)

def process_single_file(vcf_folder, genome_name, verbose=False):
    vcf_file_gz = os.path.join(vcf_folder, f"{genome_name}.h.vcf.gz")
    vcf_file_plain = os.path.join(vcf_folder, f"{genome_name}.h.vcf")
    
    if os.path.exists(vcf_file_gz):
        vcf_file = vcf_file_gz
    elif os.path.exists(vcf_file_plain):
        vcf_file = vcf_file_plain
    else:
        print(f"ERROR: VCF file not found for {genome_name}")
        return

    bed_file = os.path.join(vcf_folder, f"{genome_name}.h.bed")

    try:
        if vcf_file.endswith(".gz"):
            with gzip.open(vcf_file, "rt") as f:
                output = parse_vcf_lines(f, genome_name, verbose=verbose)
        else:
            with open(vcf_file, "r") as f:
                output = parse_vcf_lines(f, genome_name, verbose=verbose)

        # Sort output
        output_sorted = sorted(output, key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))

        with open(bed_file, "w") as out:
            out.write("#chrom\tstart\tend\tstrand\tchecksum\tgenome\tref_chr\tref_start\tref_end\tref_checksum\n")
            out.write("\n".join(output_sorted) + "\n")

        print(f"DONE {genome_name}: {len(output_sorted)} haplotype blocks")

    except Exception as e:
        print(f"ERROR processing {genome_name}: {e}")

if __name__ == "__main__":
    main()