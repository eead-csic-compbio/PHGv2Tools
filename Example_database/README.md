To guide any user about the usage of phgtools an example set has been prepared. It uses same files as [PHG testing kit]
It actually takes: 
- LineA.fa, LineB.fa, LineC.fa
- Ref.fa
- anchors.gff

*** Ensure to have the proper phgtools conda and installed the phg software

It will give a 3 genomes pangenome using a reference genome and its gene annotation.
There is also a raw fastq (doomie) imitating a low density sequence.

This is the initial material to start the pangenome:
```
phg_v2_example/
├── data
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA.fa
│   ├── LineB.fa
│   ├── LineC.fa
│   └── LineC.fastq
└── output
    └── vcf_files
```
After installing PHG, ensure to [initiate the database](https://github.com/maize-genetics/phg_v2?tab=readme-ov-file#build-and-load-data) and to [prepare properly the fastas](https://phg.maizegenetics.net/build_and_load/#prepare-assembly-fasta-files).

```
./phg initdb --db-path /path/to/dbs
./phg prepare-assemblies --keyfile /path/to/keyfile --output-dir data/updated_assemblies --threads numberThreadstoRun
```
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│    └── vcf_files
└── vcf_dbs
    └── gvcf_dataset
    └── hvcf_dataset
    └── temp
```
Then, prepare the ranges reference file, align and compress the genomes to include at the database. Follow PHG commands
```
./phg create-ranges --reference-file data/updated_assemblies/Ref.fa --gff my.gff --boundary gene --pad 500 --range-min-size 500 -o /path/to/bed/file.bed
./phg align-assemblies --gff anchors.gff --reference-file data/updated_assemblies/Ref.fa --assembly-file-list assembliesList.txt --total-threads 20 --in-parallel 4 -o /path/for/generatedFiles
./phg agc-compress --db-path /path/to/dbs --reference-file data/updated_assemblies/Ref.fa --fasta-list /my/assemblyFastaList.txt
```
Directory now should be:
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│   │   └── vcf_files
│   ├── anchorwave_gff2seq_error.log *
│   ├── anchorwave_gff2seq_output.log *
│   ├── LineA.maf *
│   ├── LineA.sam *
│   ├── LineA.svg *
│   ├── LineA_Ref.anchorspro *
│   ├── LineB.maf *
│   ├── LineB.sam *
│   ├── LineB.svg *
│   ├── LineB_Ref.anchorspro *
│   ├── LineC.maf *
│   ├── LineC.sam *
│   ├── LineC.svg *
│   ├── LineC_Ref.anchorspro *
│   ├── minimap2_LineA_error.log *
│   ├── minimap2_LineA_output.log *
│   ├── minimap2_LineB_error.log *
│   ├── minimap2_LineB_output.log *
│   ├── minimap2_LineC_error.log *
│   ├── minimap2_LineC_output.log *
│   ├── minimap2_Ref_error.log *
│   ├── minimap2_Ref_output.log *
│   ├── proali_LineA_outputAndError.log *
│   ├── proali_LineB_outputAndError.log *
│   ├── proali_LineC_outputAndError.log *
│   ├── ref.cds.fasta *
│   └── Ref.sam *
└── vcf_dbs
    ├── assemblies.agc*
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── reference *
    │   ├── Ref.bed *
    │   └── Ref.fa *
    └── temp/
```
Then, proceed creating the VCF files
```
./phg create-ref-vcf --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --reference-url https://url-for-ref --reference-name B73 --db-path /path/to/tiled/dataset folder
./phg create-maf-vcf --db-path /path/to/dbs --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --maf-dir /my/maf/files -o /path/to/vcfs
```
New files are generated.
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│   │   └── vcf_files
│   │       ├── LineA.h.vcf.gz *
│   │       ├── LineA.h.vcf.gz.csi *
│   │       ├── LineA.g.vcf.gz *
│   │       ├── LineA.g.vcf.gz.csi *
│   │       ├── LineB.h.vcf.gz *
│   │       ├── LineB.h.vcf.gz.csi *
│   │       ├── LineB.g.vcf.gz *
│   │       ├── LineB.g.vcf.gz.csi *
│   │       ├── LineC.h.vcf.gz *
│   │       ├── LineC.h.vcf.gz.csi *
│   │       ├── LineC.g.vcf.gz *
│   │       ├── LineC.g.vcf.gz.csi 
│   │       └── VCFMetrics.tsv *
│   ├── anchorwave_gff2seq_error.log 
│   ├── anchorwave_gff2seq_output.log 
│   ├── LineA.maf 
│   ├── LineA.sam 
.   .
.   .
.   .
│   └── Ref.sam 
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files *
    │   ├── Ref.h.vcf.gz *
    │   └── Ref.h.vcf.gz.csi *
    ├── reference 
    │   ├── Ref.bed 
    │   └── Ref.fa 
    └── temp/
```
Then, building the kmer index and prepare for the imputation of the test fastq file
```
./phg build-kmer-index --db-path /my/db/uri --hvcf-dir /my/hvcf/dir
./phg map-kmers --hvcf-dir /my/hvcf/dir --kmer-index /my/hvcf/dir/kmerIndex.txt --key-file /my/path/keyfile --output-dir /my/mapping/dir
./phg find-paths --path-keyfile /my/path/keyfile --hvcf-dir /my/hvcf/dir --reference-genome /my/ref/genome --path-type haploid --output-dir /my/imputed/hvcfs
```
After imputation, the new LineD.h.vcf represents the imputed haplotype.
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│   │   └── vcf_files
│   │       ├── LineA.h.vcf.gz 
│   │       ├── LineB.h.vcf.gz 
│   │       ├── LineC.h.vcf.gz 
│   │       .
│   │       .
│   │       ├── kmerIndex.txt *
│   │       ├── kmerIndexStatistics.txt *
│   │       ├── LineD_readMapping.txt *
│   │       ├── LineD_imputed.h.vcf *
│   │       └── VCFMetrics.tsv 
│   ├── anchorwave_gff2seq_error.log 
│   ├── anchorwave_gff2seq_output.log 
│   ├── LineA.maf 
│   ├── LineA.sam 
.   .
.   .
.   .
│   └── Ref.sam 
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files *
    │   ├── Ref.h.vcf.gz *
    │   └── Ref.h.vcf.gz.csi *
    ├── reference 
    │   ├── Ref.bed 
    │   └── Ref.fa 
    └── temp/
```

With samtools, index the genomes fasta files for further cropping of sequences.
The following oneliner prepare all the genomes of the pangenome fastas folder:
```
perl -e 'for my $file (glob("/scratch/PHG/PHGv2Tools/Example_database/data/prepared_genomes/*.fa")) { print "Indexing $file\n"; system("samtools faidx $file") == 0 or die "Failed to index $file: $!"; }'
```

Last step, to merge all hvcfs of Pangenome (LineA, LineB, LineC) in a single file using:
```
phg merge-hvcfs --input-dir my/hvcf/directory --output-file output/merged_hvcfs.h.vcf --id-format CHECKSUM --reference-file --range-bedfile
```
Instead of merged hvcf files we can use a range summary table for utilities as ```fasta-from-key```.
Run the phg command:
```
./phg hapid-sample-table --hvcf-dir <Path to directory holding hVCF files>  --output-file <.tsv output full path>
```

This is then the files to have in mind for phgtools analysis:

```
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── Ref-v5.fai
│   │   ├── LineA.fa
│   │   ├── LineA.fai
│   │   ├── LineB.fa
│   │   ├── LineB.fai
│   │   ├── LineC.fa
│   │   └── LineC.fai
└── output
│   │   └── vcf_files
│   │       ├── LineA.h.vcf.gz ***
│   │       ├── LineB.h.vcf.gz ***
│   │       ├── LineC.h.vcf.gz ***
│   │       ├── LineD_imputed.h.vcf ***
│   │       ├── hapIDtable.tsv
│   │       └── MergedLinesA_B_C.h.vcf ***
└── vcf_dbs
    └── hvcf_files 
            └──── Ref.h.vcf.gz ***
```
## PHG tools utilities:

#### Pangenome ranges evolution
```
range-pangenome-evolution --hvcf-folder /Example_database/output/vcf_files/ --reference-file Ref.fa --range-bedfile output/ref_ranges.bed
```
The testing files are small, and new ranges are not included after including new genomes
![RangesAmplificationSlope](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/RangesAmplificationSlope.png)


#### Core, accessory and unique ranges
```
core-range-detecter --pangenome-hvcf output/MergedLinesA_B_C.h.vcf
```
These pangenome has either ranges in all haplotypes or only in one, as shown in plots:

![MergedLinesA_B_C.2.png](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/MergedLinesA_B_C.h.vcf_2.png)
![MergedLinesA_B_C.2.png](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/MergedLinesA_B_C.h.vcf_1.png)

#### Plot pangenome region/chromosomes
If not region is specified, whole chromosome will be plotted
```
plot-pangenome-chromosomes --hvcf-folder output/vcf_files/ --reference-hvcf output/Ref.h.vcf.gz --chromosome chr1 --reference-fasta Ref.fa 
```
![chr1_FULL](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/pangenome_FULL_chr1.png)
On the other hand:
```
plot-pangenome-chromosomes --hvcf-folder output/vcf_files/ --reference-hvcf output/Ref.h.vcf.gz --chromosome chr2 --region 15000-35000 --reference-fasta Ref.fa 
```
![chr2_15000-35000](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/pangenome_chr2_15000-35000.png)

#### Check haplotypes for a region
```
check-haplotype-alleles --hvcf output/MergedLinesA_B_C.h.vcf --reference-fasta Ref.fa --start 18800 --end 20100 --chromosome 2
```
It will give as output a list of the ranges existing inside your coordinates, or in which range is contained:
```
there are 2 keys of ranges containing the coordinates:
['3efc16790e55a2a8334c939d0795dfde', '6727268d7ef78af2591fe7be98470895']
Line: LineA     Region: 2:17607-21500   Reference range:2:18001-21500> 
Line: LineB     Region: 2:18001-21500   Reference range:2:18001-21500>
```

#### Check Identity against pangenome
```
check-imputated-haplotype --hvcf-folder output/ --hvcf-file output/LineD.h.vcf
```
It gives all matches with pangenome haplotypes against the imputed one:
```
Total ranges in output/LineD.h.vcf is 38
Checking <_io.TextIOWrapper name='output/LineB.h.vcf.gz' encoding='UTF-8'>
Match count for output/LineB.h.vcf.gz is 0 out of 38 ranges(0.0%)
Checking <_io.TextIOWrapper name='output/LineA.h.vcf.gz' encoding='UTF-8'>
Match count for output/LineA.h.vcf.gz is 38 out of 38 ranges(100.0%)
Checking <_io.TextIOWrapper name='output/LineC.h.vcf.gz' encoding='UTF-8'>
Match count for output/LineC.h.vcf.gz is 21 out of 38 ranges(55.26%)
```

#### Plot imputed genome
```
plot-imputed-hvcf --input-hvcf output/LineD.h.vcf --pangenome-hvcf-folder output/ --reference-hvcf hvcf_files/Ref.h.vcf.gz
```
It will plot the imputed h.vcf file:
![LineD_hvcf_plot](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/LineD.h.vcf.png)

#### Export fasta sequence from range

Note: before starting, ensure to have the merged h.vcf file or the hapIDtable.tsv file.

```
fasta-from-key --fastas-folder data/prepared_genomes --vcf-folder output/vcf_files/ --key 3efc16790e55a2a8334c939d0795dfde --merged-vcf output/merged_hvcfs.h.vcf
```
Or use instead:
```
fasta-from-key --fastas-folder data/prepared_genomes --vcf-folder output/vcf_files/ --key 3efc16790e55a2a8334c939d0795dfde --hapIDtable output/vcf_files/hapIDtable.tsv
```
If also adding ```--output output/vcf_files/``` will result in:
```
opening output/vcf_files//LineA.h.vcf.gz to extract the coordinates

Sample name:  LineA
Chromosome:  2
Start:  21501
End:  23500

______________________________________________________________________________________

cropping sequence for LineA from 2:21501-23500

output file: output/vcf_files//LineA_2_21501_23500.fa
```
And the file will contain:
``` 
>2:21501-23500
GCACGGCACGAAAGGAGCGAACCCGAGGGCGGGGGACAAGGAAGACCGACCCCGGACAAC
CCCACCCCGGAACAGAGCGAGACGGGAACTAGAGGCGGCCCCACCGGGGGCAACCTGGAA
[...]
ATACCAGCGCGGAAGAACACGAGGCAAAAGGGACGAACCACGAGCAAAACGCCGGAGAAG
CCGCACAGAGCCAAGCGGGGCGAGCCGAACGAAACCGCAGGGCGCCAACGGGGAAACAGG
GCGACAGCAACGACGGACACGGCCCGGAAAGGCACAACCAGGAGCGCCACCCCAGGACGC
AGCCAGAGACCGGCAACCCGGGGGGGGAGGAGGCCCCGGGCGTCGCCGCACGACGACACT
CCAAGAGACCAGGCGATAGAGCGGGACGAATGAGAAGCCGCACAAACGACGAACGCGACC
AAAAAGCGGAGGGGCAAAGGGCCACCACGGGAGCAGCACCACCAGGGGACCGAAACAGGG
ACATCCAACACACGGGAGCG
```