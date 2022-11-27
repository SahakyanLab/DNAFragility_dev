#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path=$3
cores=$4
interval=$5
alignment_strands=$6

# create directories to save files
mkdir -p ../../Raw_data/$breakpoint_experiment/breakpoint_positions/chr{1..22}
mkdir -p ../../Raw_data/$breakpoint_experiment/breakpoint_positions/concat
mkdir -p ../../Raw_data/$breakpoint_experiment/average_levdist/
mkdir -p ../figures/$breakpoint_experiment

# check if any BAM files present
nr_of_BAM_files=$((ls ../../Raw_data/$breakpoint_experiment/*.bam | wc -l) 2> /dev/null)
nr_of_fasta_files=$((ls ../../Raw_data/$breakpoint_experiment/*.fasta.gz | wc -l) 2> /dev/null)

# Process files
Rscript Process.R $my_path