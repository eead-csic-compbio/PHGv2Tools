#!/usr/bin/env python
# coding: utf-8

# Script to extract the fasta sequence from a key of a PHGv2 pangenome.
# Prepared to work either from the merged h.hvcf file or HapIDtable:
# 
# -merged h.vcf: $ phg merge-hvcf
# -HapIDtable: $ phg hapid-sample-table

# In[16]:


def AmIaNotebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type of environment
    except NameError:
        return False  # Not in an interactive environment

#Check if the name of this file ends with .py or with .


# In[17]:


def ExtractCoordFromKey(merged_vcf, hapIDtable, key, vcf_folder):

    if merged_vcf is not None:
        with open(merged_vcf, "r") as f:
            for line in f:
                if line.startswith("##ALT"):
                    if key in line:
                        #print(line) #debug
                        sample_name = line.split(",")[3]
                        sample_name = sample_name.split("=")[1]
                        region = line.split(",")[4]
                        region = region.split("=")[1]
                        #print(region) #debug
                        chr = region.split(":")[0]
                        region = region.split(":")[1]
                        start = region.split("-")[0]
                        end = region.split("-")[1]
                        break

    elif hapIDtable is not None:
        with open(hapIDtable, "r") as f:
            for line in f:
                if key in line:
                    #print(line) #debug
                    sample_name = line.split("\t")[1]
                    sample_name = sample_name.split(",")[0] #only keep the first genome with the key, enough
                    hvcf_file = f"{vcf_folder}/{sample_name}.h.vcf.gz"
                    print(f"opening {hvcf_file} to extract the coordinates") #debug
                    with gzip.open(hvcf_file, "rt") as f:
                        for line in f:
                            if line.startswith("##ALT"):
                                continue
                            elif key in line:
                                #print(line) #debug
                                chr = line.split("\t")[0]
                                start = line.split("\t")[1]
                                end = line.split("\t")[7]
                                end = end.split("=")[1]
                    break

    print("Sample name: ", sample_name)
    print("Chromosome: ", chr)
    print("Start: ", start)
    print("End: ", end)

    return sample_name, chr, start, end


# In[18]:


def ExtractFastaFromCoord(fastas_folder, sample_name, chr, start, end, output_folder):
    fasta_file = os.path.join(fastas_folder, sample_name + ".fa")

    print("\n______________________________________________________________________________________\n")

    print (f"cropping sequence for {sample_name} from {chr}:{start}-{end}\n")

    #try to open the fasta file
    try:
        with open(fasta_file, "r") as f:
            pass
    except FileNotFoundError:
        print("The fasta file for the sample could not be found.")
        exit()

    #Use samtools faidx to extract the sequence
    os.system(f"samtools faidx {fasta_file} {chr}:{start}-{end} > {sample_name}_{chr}_{start}_{end}.fa")

    #print the sequence

    if output_folder is None:
        #print it at terminal
        with open(f"{sample_name}_{chr}_{start}_{end}.fa", "r") as f:
            for line in f:
                print(line)

        os.system(f"rm {sample_name}_{chr}_{start}_{end}.fa")

    elif output_folder is not None:
        #print it in a file
        with open(f"{output_folder}/{sample_name}_{chr}_{start}_{end}.fa", "w") as out, open (f"{sample_name}_{chr}_{start}_{end}.fa", "r") as f:
            for line in f:
                out.write(line)

        print(f"output file: {output_folder}/{sample_name}_{chr}_{start}_{end}.fa")


# In[ ]:


import os
import argparse
import gzip

if AmIaNotebook() == True:
# Edit here the stuff if you are working here as a notebook
    vcf_folder = "/scratch/PHG_panbarleyV1/vcf_dbs/hvcf_files/"
    merged_vcf = "/scratch/PHG_panbarleyV1/output/vcf_files/merged_08112024.h.vcf"
    fastas_folder = "/scratch/PHG_panbarleyV1/data/prepared_assemblies/"
    hapIDtable = "/scratch/PHG_panbarleyV1_test/output/graph/hapIDtable_121224.tsv"
    output_folder = "/scratch/PHG_panbarleyV1_test/output/graph/" #If None, output will be just printed at the console
    key = input("Please input the key to extract the fasta:")

def main(merged_vcf, fastas_folder, vcf_folder, hapIDtable, key, output_folder):

    """
    This script is used to extract fasta sequence from a range key of a pangenome build with PHGv2
    It can use either a merged vcf file or a haplotype ID table to extract the fasta
    Requires the path of the folder containing the vcf files, the path of the folder containing the fasta files, the key to extract the fasta,
    end either the merged hvcf file or the haplotype ID table, builts with merge-hvcf or hapid-sample-table options of the PHGv2 pipeline
    If the output folder is not specified, the fasta will be printed at the console
    If there are more than one genomes in pangenome owning the key, the fasta will be extracted from the first genome found (It doesn't matter, all are the same sequences)
    
    !NOTE: It uses an indexed fasta with samtools faidx. If the fasta is not indexed, it will not work.
    Check annotation of PHGtools for more information about the fasta indexing.
    
    It can be used in a notebook or in a terminal
    """

    if merged_vcf is None and hapIDtable is None:
        print("Please input the merged vcf file or haplotype ID table. At least one of them")
        exit(1) 
        # merged_vcf or hapIDtable can be one of them None, but not both
    
    sample_name, chr, start, end = ExtractCoordFromKey(merged_vcf, hapIDtable, key, vcf_folder)
    ExtractFastaFromCoord(fastas_folder, sample_name, chr, start, end, output_folder)

    



# In[20]:


if __name__ == "__main__":
    try:
        main(vcf_folder, merged_vcf, fastas_folder, hapIDtable, output_folder, key)
    except KeyboardInterrupt:
        print("Script interrupted by user. Exiting...")
        sys.exit(0)
        raise KeyboardInterrupt

