#!/usr/bin/env python
# coding: utf-8

# # Code to process a map_kmers file output of phgv2. It plots both the % of identity of each genome to the query genome but also plots it with pyGenomeTracks
# 
# # Perfomed from output of 
# # phg version 2.4.4.158
# 
# # Code developed at 31/04/2025
# 
# # Run with phgtools conda env, python, using 
# 
# # Uses as imput the map_kmers file, and the hapIDranges file output of phg sample-hapid-by-range

# In[ ]:


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

#Check if the name of this file ends with .py or with .ipynb


# In[ ]:


def DetermineGenomeSample(hapID_ranges_file, query_genome):

    with open (hapID_ranges_file, 'r') as f:
        header = f.readline().strip()
        header = header.split('\t')
        header_hash = {}
        for i in range(len(header)):
            header_hash[header[i]] = i

    print ("Header hash: ", header_hash)

    #check the column of the query genome
    if query_genome not in header_hash:
        print ("Error: ", query_genome, " not found in the hapID ranges file")
        return None
    else:
        query_genome_index = header_hash[query_genome]
        print (f"Genome {query_genome} is the one at the column: ", query_genome_index, "in the hapID ranges file")

    return query_genome_index




# In[ ]:


def GenerateTempBed (hapID_ranges_file, query_genome_index):

    # Get the folder of the map kmers file
    hapID_ranges_file_folder = hapID_ranges_file.split("/")[:-1]
    hapID_ranges_file_folder = "/".join(hapID_ranges_file_folder)
    hapID_ranges_file_folder = hapID_ranges_file_folder + "/"

    #generate a temporary .bed file
    hapID_ranges_file_bed = hapID_ranges_file_folder + "temp.bed"

    with open (hapID_ranges_file, 'r') as f, open (hapID_ranges_file_bed, 'w') as f_bed:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line == "":
                continue
            line = line.split("\t")
            if line[query_genome_index] == ".":
                continue
            else:
                chr = line[0]
                start = line[1]
                end = line[2]
                key = line[query_genome_index]
            # Write to the bed file
                f_bed.write (f"{chr}\t{start}\t{end}\t{key}\n")

    return hapID_ranges_file_bed



# In[ ]:


def CheckIntersectionBed_Mapping(hapID_ranges_file_bed, map_kmers_file_bed):   
    from tqdm import tqdm

    # Generate a temp_2.bed file
    hapID_ranges_file_bed_query = hapID_ranges_file_bed.split("/")[:-1]
    hapID_ranges_file_bed_query = "/".join(hapID_ranges_file_bed_query)
    hapID_ranges_file_bed_query = hapID_ranges_file_bed_query + "/temp_2.bed"

    # Load map_kmers_file_bed into memory as a dictionary for fast lookups
    # If not, It used to took hours lol
    map_kmers_dict = {}
    with open(map_kmers_file_bed, 'r') as map_kmers:
        for line in map_kmers:
            if line.startswith("#") or line.strip() == "" or line.startswith("HapIds"):
                continue
            line = line.strip().split("\t")
            map_keys = line[0].split(",")  # Extract keys
            for key in map_keys:
                map_kmers_dict[key] = line  # Store the entire line for reference

    # Count the total number of lines in hapID_ranges_file_bed for the progress bar
    with open(hapID_ranges_file_bed, 'r') as hapID:
        total_lines = sum(1 for _ in hapID)

    # Process hapID_ranges_file_bed and write matches to temp_2.bed
    with open(hapID_ranges_file_bed, 'r') as hapID, open(hapID_ranges_file_bed_query, 'w') as final_bed:
        for line in tqdm(hapID, total=total_lines, desc="Processing hapID ranges"):
            line = line.strip().split("\t")
            if len(line) < 4:  # Skip lines that don't have enough columns
                continue
            key = line[3]
            # Remove the <> of the key
            key = key.replace("<", "").replace(">", "")

            # Check if the key exists in the map_kmers_dict
            if key in map_kmers_dict:
                final_bed.write(f"{line[0]}\t{line[1]}\t{line[2]}\t{key}\n")

    return hapID_ranges_file_bed_query


# In[ ]:


def CaptureChrFromBED (hapID_ranges_file_bed_query):
    """
    This function will capture the chromosome from the BED file, storing them in a list
    """
    chr_list = []
    with open (hapID_ranges_file_bed_query, 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
        for line in lines:
            chr = line.split('\t')[0]
            if chr not in chr_list:
                chr_list.append(chr)
    chr_list.sort()
    #print(f"Chromosomes found in the BED file are: {chr_list}\n")
    return chr_list


# In[ ]:


def BEDSplitByChr (hapID_ranges_file_bed_query):
    pass
    """
    splits and sorts the BED file in different files by chromosome, sorting it by start position. 
    This is only in order to plot it with pygenometracks. Takes by output the BED file generated by hvcf2bed.

    """

    chr_list = CaptureChrFromBED(hapID_ranges_file_bed_query)
    print(f"Chromosomes found in the BED file are: {chr_list}\n")
    if len(chr_list) == 0:
        print(f"Chromosome list is empty, exiting")
        exit()

    plot_folder = hapID_ranges_file_bed_query.split("/")[:-1]
    plot_folder = "/".join(plot_folder)
    plot_folder = plot_folder + "/"
    #print("Command CaptureChrFromBED executed")

    for chr in chr_list:    #check if the files splited by chr already exist, and if they do, remove them            
        chr_file = (f"{plot_folder}TEMP_{os.path.basename(hapID_ranges_file_bed_query).split('.')[0]}_{chr}.bed")
        if os.path.exists(chr_file):
            print(f"File {chr_file} already exists, removing it and building it again") 
            os.remove (chr_file)

        try:
            open (hapID_ranges_file_bed_query, 'r')
        except FileNotFoundError:
            print(f"File {hapID_ranges_file_bed_query} not found")

        with open (hapID_ranges_file_bed_query, 'r') as f:
            lines = f.readlines()
            if len(lines) == 0:
                raise Exception(f"File {hapID_ranges_file_bed_query} is empty")
            if len(chr_list) == 0:
                raise Exception(f"Chromosome list is empty")

            #print(f"opening {hvcf_BED} to split it by chromosome")

            for chr in chr_list:
                chr_file = (f"{plot_folder}TEMP_{os.path.basename(hapID_ranges_file_bed_query).split('.')[0]}_{chr}.bed")
                with open (chr_file, 'w') as chr_bed:
                    chr_bed.write("chrom\tchromStart\tchromEnd\tname\n")
                    for line in lines:
                        if chr in line:
                            #substitute {chr} with Chr
                            line = line.replace(chr, f"Reference_genome")
                            chr_bed.write(line)
                try:
                    open(chr_file, 'r')
                except FileNotFoundError:
                    print(f"File {chr_file} not found") 
                else:
                    if len(open(chr_file).readlines()) <= 1:  
                        raise Exception(f"File {chr_file} of the {chr} is empty\n Im actually trying to find the chromosome {chr} in the file {hvcf_BED}")
            #        print(f"Finished splitting {hvcf_BED} by chromosome to {chr}")

    return chr_list


# In[ ]:


def PrepareTrackPlot (query_genome, hapID_ranges_file_bed_query, chr_list):

    plot_folder = hapID_ranges_file_bed_query.split("/")[:-1]
    plot_folder = "/".join(plot_folder)
    plot_folder = plot_folder + "/"


    track_ini = (f"{plot_folder}TEMP_{query_genome}_track.ini")

    if os.path.exists(track_ini):
        print(f"File {track_ini} already exists, removing it and building it again") 
        os.remove(track_ini)

    with open (track_ini, 'w') as track:

        for chr in chr_list:
            track.write(f"[haplotype_{chr}]\n")
            track.write(f"title = {chr}\n")
            track.write(f"file = {plot_folder}TEMP_temp_2_{chr}.bed \n")
            track.write(f"height = 3\n")
            track.write(f"color = bed_rgb\n")
            track.write("display = collapsed\n")
            track.write("labels = false\n")
            # track.write("label_fontsize = 10\n") #unused
            track.write("border_color = black\n")
            track.write("line_width = 0.1\n")  # Changed border size
            track.write("\n")
            track.write("[spacer]\nheight = 0.5\n")

        track.write("[x-axis]\n")
        track.write("where = bottom\n")
        #track.write("label = true\n")
        #track.write("font_size = 10\n")
        #track.write(f"title = {query_genome} using as reference MorexV3\n")
        track.write("[spacer]\nheight = 0.5\n")

    print(f"Finished preparing the tracks for pygenometracks at {track_ini}\n")
    return track_ini


# In[ ]:


def SetExtremesToPlot(hapID_ranges_file_bed_query, chr_list):
    #capturing the extreme values to plot of each chromosome

    plot_folder = hapID_ranges_file_bed_query.split("/")[:-1]
    plot_folder = "/".join(plot_folder)
    plot_folder = plot_folder + "/"

    chr_max = None
    for chr in chr_list:
        chr_file = (f"{plot_folder}TEMP_temp_2_{chr}.bed")
        with open (chr_file, 'r') as f:
            #print (f"Processing {chr_file}")
            lines = f.readlines()
            if lines == []:
                raise Exception(f"Empty file {chr_file}")

            try:    
                last_line = lines[-1]
            except:
                print(f"This is the file im working now: {lines}\n")
                raise IndexError(f"I can not find the last line in {chr_file}")

            last_line = lines[-1]
            last_coord = last_line.split('\t')[2]
            if chr_max == None or chr_max < last_coord:
                chr_max = last_coord
    chr_max = int(int(chr_max)*1.005)
    print(f"Max coord to plot will be {chr_max}\n")

    return chr_max


# In[ ]:


def RunPygenometracksCommandHVCF(chr_max, track_ini, chr_list, query_genome, hapID_ranges_file_bed_query):

    """
    Building and executing the command for plotting:
    It calls the functions SetExtremesToPlot and PrepareTrack_HvcfPlot
    """

    plot_folder = track_ini.split("/")[:-1]
    plot_folder = "/".join(plot_folder)
    plot_folder = plot_folder + "/"

    height = len(chr_list)*10

    command = (f"pyGenomeTracks --tracks {track_ini} --region Reference_genome:1-{chr_max} --outFileName {plot_folder}/{os.path.basename(hapID_ranges_file_bed_query)}.png --dpi 500")

    try: 
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"running the command: {command}\n")
        #print(result.stdout)
        #print(result.stderr)
        if "No valid intervals were found" in result.stdout or "No valid intervals were found" in result.stderr:
            raise Exception(f"Command {command} failed with error: {result.stderr}")
    except subprocess.CalledProcessError as e:
        raise Exception(f"Command {command} failed with error: {e.stderr}")



# In[ ]:


def CalculatePresentRanges(hapID_ranges_file_bed, hapID_ranges_file_bed_query):

    #Divide the amount of lines in the hapID_ranges_file_bed by the amount of lines in the hapID_ranges_file_bed_query

    with open (hapID_ranges_file_bed, 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
        hapID_ranges_file_bed_lines = len(lines)
        #print(f"Lines in {hapID_ranges_file_bed}: {hapID_ranges_file_bed_lines}")
    with open (hapID_ranges_file_bed_query, 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
        hapID_ranges_file_bed_query_lines = len(lines)
        #print(f"Lines in {hapID_ranges_file_bed_query}: {hapID_ranges_file_bed_query_lines}")
    if hapID_ranges_file_bed_lines == 0:
        raise Exception(f"File {hapID_ranges_file_bed} is empty")
    if hapID_ranges_file_bed_query_lines == 0:
        raise Exception(f"File {hapID_ranges_file_bed_query} is empty")
    if hapID_ranges_file_bed_lines == hapID_ranges_file_bed_query_lines:
        print(f"All the ranges are present in the hapID_ranges_file_bed_query")
    else:
        present_percentage = (hapID_ranges_file_bed_query_lines / hapID_ranges_file_bed_lines) * 100
        print(f"Present ranges: {hapID_ranges_file_bed_query_lines} / {hapID_ranges_file_bed_lines} ({present_percentage:.2f}%)")
    #return hapID_ranges_file_bed_lines / hapID_ranges_file_bed_query_lines


# In[ ]:


import os
import gzip
import re
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
from matplotlib.lines import Line2D
import argparse
import sys
from tqdm import tqdm

if AmIaNotebook() == True:
    map_kmers_file = "/scratch/PHG_v2/output/ERR753323_readMapping.txt"
    hapID_ranges_file = "/scratch/PHG_v2/output/hapIDranges.tsv"
    query_genome = "MorexV3"
    skip_plot = False

def main(map_kmers_file, hapID_ranges_file, query_genome, skip_plot):

    """
        Main function to run the script: checks, for a single genome of the pangenome, how many ranges have been found in the mapping of a fastq file to the pangenome
        It uses the map_kmers file, and the hapIDranges file output of phg sample-hapid-by-range
        It generates a plot with pygenometracks
    """

    print(f"Running the script with the following parameters:\n")
    print(f"Map kmers file: {map_kmers_file}\n")
    print(f"HapID ranges file: {hapID_ranges_file}\n")
    print(f"Query genome: {query_genome}\n")    
    print(f"Skip plot: {skip_plot}\n")


    if not os.path.exists(map_kmers_file):
        print(f"File {map_kmers_file} does not exist")
        exit()
    if not os.path.exists(hapID_ranges_file):
        print(f"File {hapID_ranges_file} does not exist")
        exit()
    if query_genome == None:
        print(f"Please provide a reference genome")
        exit()
    if query_genome == "":
        print(f"Please provide a reference genome")
        exit()

    query_genome_index = DetermineGenomeSample(hapID_ranges_file, query_genome)
    hapID_ranges_file_bed = GenerateTempBed(hapID_ranges_file, query_genome_index)
    hapID_ranges_file_bed_query = CheckIntersectionBed_Mapping(hapID_ranges_file_bed, map_kmers_file)

    if not skip_plot:
        chr_list = BEDSplitByChr(hapID_ranges_file_bed_query)
        track_ini = PrepareTrackPlot(query_genome, hapID_ranges_file_bed_query, chr_list)
        chr_max = SetExtremesToPlot(hapID_ranges_file_bed_query, chr_list)
        RunPygenometracksCommandHVCF(chr_max, track_ini, chr_list, query_genome, hapID_ranges_file_bed_query)

    CalculatePresentRanges(hapID_ranges_file_bed, hapID_ranges_file_bed_query)


    #remove temp files
    plot_folder = hapID_ranges_file_bed_query.split("/")[:-1]
    plot_folder = "/".join(plot_folder)
    plot_folder = plot_folder + "/"

    os.remove(hapID_ranges_file_bed)
    os.remove(hapID_ranges_file_bed_query)
    if skip_plot == False:
        os.remove(track_ini)
        for chr in chr_list:
            chr_file = (f"{plot_folder}TEMP_temp_2_{chr}.bed")
            os.remove(chr_file)
    print(f"\nFinished removing temp files\n")

    #rename the output file to query_name_and the map_kmers_file
    if skip_plot == False:
        output_png = f"{plot_folder}{os.path.basename(hapID_ranges_file_bed_query)}.png"
        os.rename(output_png, f"{plot_folder}{query_genome}_{os.path.basename(map_kmers_file)}_{query_genome}.png")
        print (f"The output file is: {plot_folder}{query_genome}_{os.path.basename(map_kmers_file)}_{query_genome}.png\n")


# In[ ]:


if __name__ == "__main__":

    main(map_kmers_file, hapID_ranges_file, query_genome, skip_plot)

