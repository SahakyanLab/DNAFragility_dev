#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path=$3
cores=$4
chromosome=$5
control=$6
category=$7
auto_fit=$8

# create directory to save RMSD results
mkdir -p ../{data,figures}/${breakpoint_experiment}

for kmer in 2 4 6 8 10
do
	echo "Calculating RMSD for kmer $kmer..."
	Rscript RMSD.R $my_path $breakpoint_experiment $chromosome $ref_path $kmer $cores $control
done

echo "Plotting RMSD figures..."
Rscript RMSDPlots.R $my_path $breakpoint_experiment $chromosome $category $auto_fit