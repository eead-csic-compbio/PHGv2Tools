import argparse
from PHGv2Tools.Modules.CheckHaplotypeAllelesInPangenome import main as check_haplotype_alleles_main
from PHGv2Tools.Modules.CheckImputatedHaplotype import main as check_imputated_haplotype_main
from PHGv2Tools.Modules.CoreRangeDetecter import main as core_range_detecter_main
from PHGv2Tools.Modules.PlotImputedHvcf import main as plot_imputed_hvcf_main
from PHGv2Tools.Modules.PlotPangenomeChromosomes import main as plot_pangenome_chromosomes_main
from PHGv2Tools.Modules.RangePangenomeEvolution import main as range_pangenome_evolution_main
from PHGv2Tools.Modules.FastaFromKey import main as fasta_from_key_main
from PHGv2Tools.Modules.CheckSetup import main as check_setup
from PHGv2Tools.__version__ import __version__

def main():
    parser = argparse.ArgumentParser(prog='phgtools', description='PHG Tools Command Line Interface. CLI for the PHGv2Tools package')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    subparsers = parser.add_subparsers(dest='command')

    # Check-haplotype-alleles
    parser_check_haplotype_alleles = subparsers.add_parser('check-haplotype-alleles', help='Check haplotype alleles in pangenome')
    parser_check_haplotype_alleles.add_argument('--pangenome-hvcf', "-phvcf", type=str, help='The hvcf file', required=True)
    parser_check_haplotype_alleles.add_argument('--reference-fasta', "-ref-fa", type=str, help='The reference fasta file', required=True)
    parser_check_haplotype_alleles.add_argument('--start', "-s", type=int, help='The start coordinate')
    parser_check_haplotype_alleles.add_argument('--end', "-e", type=int, help='The end coordinate')
    parser_check_haplotype_alleles.add_argument('--chromosome', "-chr", type=str, help='The chromosome to check')

    # Check-imputated-haplotype
    parser_check_imputated_haplotype = subparsers.add_parser('check-imputated-haplotype', help='Check imputated haplotype')
    parser_check_imputated_haplotype.add_argument("--pangenome-folder", "-pf", help="Folder with the pangenome genomes", required=True)
    parser_check_imputated_haplotype.add_argument("--imputed-hvcf", "-ihvcf", help="File with the ranges to check", required=True)

    # Core-range-detecter
    parser_core_range_detecter = subparsers.add_parser('core-range-detecter', help='Core range detecter')
    parser_core_range_detecter.add_argument('--pangenome-hvcf', '-phvcf', help='Input hvcf file', required=True)

    # Plot-imputed-hvcf
    parser_plot_imputed_hvcf = subparsers.add_parser('plot-imputed-hvcf', help='Plot imputed HVCF')
    parser_plot_imputed_hvcf.add_argument("--imputed-hvcf", "-ihvcf", help="Path to the hvcf file", required=True)
    parser_plot_imputed_hvcf.add_argument("--pangenome-hvcf-folder", "-pfolder", help="Path to the folder where the plots will be saved", required=True)
    parser_plot_imputed_hvcf.add_argument("--reference-hvcf", "-refhvcf", help="Path to the reference hvcf file", required=True)

    # Plot-pangenome-chromosomes
    parser_plot_pangenome_chromosomes = subparsers.add_parser('plot-pangenome-chromosomes', help='Plot pangenome chromosomes')
    parser_plot_pangenome_chromosomes.add_argument('--pangenome-hvcf-folder', "-pfolder", help='Folder with the haplotype VCF files', required=True)
    parser_plot_pangenome_chromosomes.add_argument('--reference-hvcf', "-refhvcf", help='Reference genome hvcf file', required=True)
    parser_plot_pangenome_chromosomes.add_argument('--chromosome', "-chr", help='Chromosome to plot. Enter: chrX (chr1, chr7...)', required=True)
    parser_plot_pangenome_chromosomes.add_argument('--region', "-reg", help='Region to plot. Add it with a "-" dividing the start-end (100000-245000). If no value is provided, by default whole chr will be plotted', required=False)
    parser_plot_pangenome_chromosomes.add_argument('--reference-fasta', "-ref-fa", help='Reference genome fasta file', required=True)

    # Range-pangenome-evolution
    parser_range_pangenome_evolution = subparsers.add_parser('range-pangenome-evolution', help='Range pangenome evolution')
    parser_range_pangenome_evolution.add_argument('--pangenome-hvcf-folder', "-pfolder", type=str, help='The folder with the hvcf files', required=True)
    parser_range_pangenome_evolution.add_argument('--reference-fasta', "-ref-fa", type=str, help='The reference file', required=True)
    parser_range_pangenome_evolution.add_argument('--range-bedfile', "-bed", type=str, help='The range bedfile', required=True)

    #Fasta-from-key
    parser_fasta_from_key = subparsers.add_parser('fasta-from-key', help='Fasta from key')
    parser_fasta_from_key.add_argument('--merged-vcf', type=str, required=False, help='The path of the merged vcf file.')
    parser_fasta_from_key.add_argument('--fastas-folder', type=str, required=True, help='The path of the folder containing fasta files of pangenome.')
    parser_fasta_from_key.add_argument('--vcf-folder', type=str, required=True, help='The path of the folder containing vcf files.')
    parser_fasta_from_key.add_argument('--hapIDtable', type=str, required=False, help='The path of the haplotype ID table.')
    parser_fasta_from_key.add_argument('--key', type=str, required=True, help='The key to extract the fasta.')
    parser_fasta_from_key.add_argument('--output-folder', '-o', type=str, help='The path of the output folder.', required=False)

    # Check-setup
    parser_check_setup = subparsers.add_parser('check-setup', help='Check if phgtools, phgv2 and its dependences are callable from everywhere')

    args = parser.parse_args()


    if args.command == 'check-haplotype-alleles':
        start=None
        end=None
        chromosome=None
        check_haplotype_alleles_main(pangenome_hvcf=args.pangenome_hvcf, reference_fasta=args.reference_fasta, start=args.start, end=args.end, chromosome=args.chromosome)
    
    elif args.command == 'check-imputated-haplotype':
        check_imputated_haplotype_main(pangenome_folder=args.pangenome_folder, imputed_hvcf=args.imputed_hvcf)

    elif args.command == 'core-range-detecter':
        core_range_detecter_main(pangenome_hvcf=args.pangenome_hvcf)

    elif args.command == 'plot-imputed-hvcf':
        plot_imputed_hvcf_main(imputed_hvcf=args.imputed_hvcf, pangenome_hvcf_folder=args.pangenome_hvcf_folder, reference_hvcf=args.reference_hvcf)

    elif args.command == 'plot-pangenome-chromosomes':
        if not args.region:
            args.region = None
        plot_pangenome_chromosomes_main(pangenome_hvcf_folder=args.pangenome_hvcf_folder, reference_hvcf=args.reference_hvcf, chromosome=args.chromosome, region=args.region, reference_fasta=args.reference_fasta)
    
    elif args.command == 'range-pangenome-evolution':
        range_pangenome_evolution_main(pangenome_hvcf_folder=args.pangenome_hvcf_folder, reference_fasta=args.reference_fasta, range_bedfile=args.range_bedfile)
    
    elif args.command == 'fasta-from-key':
        # Ensure at least one of merged_vcf or hapIDtable is provided
        if not args.merged_vcf and not args.hapIDtable:
            parser.error("At least one of --merged-vcf or --hapIDtable must be provided.")
        fasta_from_key_main(merged_vcf=args.merged_vcf, fastas_folder=args.fastas_folder, vcf_folder=args.vcf_folder, hapIDtable=args.hapIDtable, key=args.key, output_folder=args.output_folder)
    
    elif args.command == 'check-setup':
        check_setup()

    else:
        parser.print_help()

if __name__ == '__main__':
    main()