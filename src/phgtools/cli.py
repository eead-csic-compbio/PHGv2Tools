import argparse
import sys
import os
import subprocess
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich import box

from phgtools import __version__

# --- IMPORT YOUR MODULES HERE ---
from phgtools.modules import hvcf2bed
from phgtools.modules import haplopainting
from phgtools.modules import fastaFromKey
from phgtools.modules import rangePangenomeEvolution
from phgtools.modules import genomeIntersection
from phgtools.modules import coreRangeDetector
from phgtools.modules import CheckHaplotypeAllelesInPangenome
from phgtools.modules import checkImputatedHaplotype
from phgtools.modules import plotImputedHvcf
from phgtools.modules import vcfDistance
from phgtools.modules import checkSetup

# --- CONFIGURATION: EDIT THIS TO CHANGE THE HELP MENU ---
COMMAND_MAP = {
    "hvcf2bed": {
        "module": hvcf2bed,
        "help": "Convert h.VCF files to BED format for visualization.",
        "usage": "phgtools hvcf2bed <folder> <genome>"
    },
    "haplopainting": {
        "module": haplopainting,
        "usage": "phgtools haplopainting --help",
        "help": "Generate haplotype painting visualizations from h.VCF files.",
    },
    "fasta-from-key": {
        "module": fastaFromKey,
        "help": "Extract FASTA sequence from a PHGv2 key/hash.",
        "usage": "phgtools fasta-from-key --key <hash> --fastas-folder <path> [--vcf-file <path> | --vcf-folder <path>] [-v]"
    },
    "range-pangenome-evolution": {
        "module": rangePangenomeEvolution,
        "help": "Analyze pangenome evolution: cumulative growth of ranges and unique haplotype keys.",
        "usage": "phgtools range-pangenome-evolution <hapIDranges.tsv> [-v]"
    },
    "genome-intersection": {
        "module": genomeIntersection,
        "help": "Analyze genome intersection from map_kmers output and BED file with haplotype keys.",
        "usage": "phgtools genome-intersection <map_kmers_file> --bed-file <bed_file> --genome-fasta <fasta_file> [-v]"
    },
    "core-range-detector": {
        "module": coreRangeDetector,
        "help": "Detect and analyze core, unique, and accessory ranges from pangenome hVCF file.",
        "usage": "phgtools core-range-detector <pangenome.hvcf> [-v] [-f <figure_format: png, pdf, svg>]"
    },
    "check-haplotype-alleles": {
        "module": CheckHaplotypeAllelesInPangenome,
        "help": "Query a pangenome hapIDranges.tsv file for overlapping genomic ranges.",
        "usage": "phgtools check-haplotype-alleles <hapIDranges.tsv> -c <chromosome> -s <start> -e <end> [-v]"
    },
    "check-imputated-haplotype": {
        "module": checkImputatedHaplotype,
        "help": "Check genomic content contribution of source genomes to an imputed haplotype (by base pair).",
        "usage": "phgtools check-imputated-haplotype <pangenome_folder> <imputed_bed_or_hvcf> [-o <output_dir>] [-v]"
    },
    "plot-imputed-hvcf": {
        "module": plotImputedHvcf,
        "help": "Plot imputed hVCF files showing genome-colored haplotype ranges across all chromosomes.",
        "usage": "phgtools plot-imputed-hvcf <pangenome_folder> <imputed_hvcf> <reference_hvcf> [-o <output_dir>] [-v]"
    },
    "vcf-distance": {
        "module": vcfDistance,
        "help": "Generate distance matrix comparing all varieties in a g.VCF file with optional heatmap visualization.",
        "usage": "phgtools vcf-distance <input.vcf.gz> [-o <out_matrix.tsv>] [-p [heatmap.pdf]] [-t <threads>] [-v]"
    },
    #
}

def print_beautiful_help():
    """Renders the Rich UI when -h or no args are passed."""
    console = Console()
    
    # 1. Title Panel
    title = Text("\n🧬 PHGv2 Command Line Tools", style="bold magenta", justify="center")
    subtitle = Text("A suite of tools for pangenome analysis \n(from PHGv2)", style="italic white", justify="center")
    
    console.print(Panel(
        Text.assemble(title, "\n", subtitle),
        box=box.ROUNDED,
        border_style="green",
        expand=False,
        padding=(1, 2)
    ), justify="center")
    console.print("") # Spacer

    # 2. Command Table
    table = Table(box=box.SIMPLE, show_header=True, header_style="bold cyan")
    table.add_column("Command", style="bold yellow", no_wrap=False)
    table.add_column("Description", style="white", no_wrap=False)
    table.add_column("Usage Hint", style="dim white", no_wrap=False)

    for cmd, info in COMMAND_MAP.items():
        table.add_row(cmd, info["help"], info["usage"])

    console.print(table)
    
    # 3. Footer
    console.print(f"\n[dim]Version: {__version__}[/dim]")
    console.print("\n[bold green]Tip:[/bold green] Run [yellow]phgtools <command> --help[/yellow] for detailed arguments.")
    console.print("[bold green]Install:[/bold green] Run [yellow]phgtools --conda-setup[/yellow] to create a conda environment with all dependencies.")
    console.print("[bold green]Setup:[/bold green] Run [yellow]phgtools --check-setup[/yellow] to verify all dependencies are installed.\n")



def run_conda_setup():
    """Run the CondaSetup.sh script to create the conda environment."""
    console = Console()
    
    # Find the CondaSetup.sh script relative to this file
    cli_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up to package root: src/phgtools -> src -> PHGv2Tools
    package_root = os.path.dirname(os.path.dirname(cli_dir))
    setup_script = os.path.join(package_root, "Misc", "CondaSetup.sh")
    
    if not os.path.exists(setup_script):
        console.print(f"[red]Error:[/red] CondaSetup.sh not found at {setup_script}")
        console.print("[dim]Make sure you have the full PHGv2Tools repository.[/dim]")
        sys.exit(1)
    
    console.print(f"[green]Running conda setup script:[/green] {setup_script}\n")
    
    # Run the script
    try:
        subprocess.run(["bash", setup_script], check=True)
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Error running setup script:[/red] {e}")
        sys.exit(1)
    except FileNotFoundError:
        console.print("[red]Error:[/red] bash not found. Please run the script manually:")
        console.print(f"[yellow]bash {setup_script}[/yellow]")
        sys.exit(1)


def main():

    # --- 1. INTERCEPT SPECIAL FLAGS ---
    # If no arguments, or explicitly asking for help, show the beautiful menu
    if len(sys.argv) == 1 or sys.argv[1] in ["-h", "--help"]:
        print_beautiful_help()
        sys.exit(0)
    
    # Handle --version
    if sys.argv[1] in ["-v", "--version"]:
        print(f"PHGv2Tools {__version__}")
        sys.exit(0)
    
    # Handle --check-setup or check-setup as a special command
    if sys.argv[1] in ["--check-setup", "check-setup"]:
        checkSetup.main(sys.argv[2:])
        sys.exit(0)
    
    # Handle conda-setup command
    if sys.argv[1] == "--conda-setup":
        run_conda_setup()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        prog="phgtools",
        description="Package to downstream analysis of pangenomes databases, working with Practical Haplotype Graph v2 and its h.VCF files"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    subparsers.required = True

# --- 2. STANDARD DISPATCHER ---
    # Disable the default help so it doesn't conflict
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest="command")

# Register commands dynamically from the map
    for cmd_name, info in COMMAND_MAP.items():
        subparsers.add_parser(
            cmd_name, 
            help=info["help"],
            add_help=False  # Important: Let the module handle its own specific --help
        )

    # Parse command (ignoring the rest of the args for now)
    args, unknown_args = parser.parse_known_args()

    # --- 3. RUN THE MODULE ---
    if args.command in COMMAND_MAP:
        try:
            # Call the main function of the selected module
            module = COMMAND_MAP[args.command]["module"]
            module.main(unknown_args)
        except AttributeError:
            print(f"Error: The module '{args.command}' does not have a 'main()' function.")

if __name__ == "__main__":
    main()
