from setuptools import setup, find_packages
import os

# Ensure the README.md file exists
if os.path.exists('README.md'):
    with open('README.md', 'r') as fh:
        long_description = fh.read()
else:
    long_description = ''

setup(
    name='PHGv2Tools',
    version='1.2.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [

            "phgtools = PHGv2Tools.phgtools:main",          # Main entry point

            # Subcommands which are actually functions in the Modules

            #"phgtools check-imputated-haplotype = Modules.CheckImputatedHaplotype:main",
            #"phgtools check-haplotype-alleles = Modules.CheckHaplotypeAllelesInPangenome:main",
            #"phgtools core-range-detecter = Modules.CoreRangeDetecter:main",
            #"phgtools plot-imputed-hvcf = Modules.PlotImputedHvcf:main",
            #"phgtools plot-pangenome-chromosomes = Modules.PlotPangenomeChromosomes:main",
            #"phgtools range-pangenome-evolution = Modules.RangePangenomeEvolution:main",
        ],
    },
    author='Joan Sarria',
    author_email='jsarria@eead.csic.es',
    description='Package to downstream analysis of pangenomes databases, working with Practical Haplotype Graph and its h.VCF files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jsarriaa/PHGv2Tools',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)