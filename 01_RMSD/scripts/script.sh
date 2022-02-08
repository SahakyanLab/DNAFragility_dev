#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path=$3
cores=$4
exp_num=$5
chromosome=$6
control=$7

# create directory to save RMSD results
mkdir -p ../{data,figures}/${breakpoint_experiment}_${exp_num}

for kmer in 2 4 6 8 10
do
	echo "Calculating RMSD for kmer $kmer..."
	Rscript RMSD.R $my_path $breakpoint_experiment $exp_num $chromosome $ref_path $kmer $cores $control
done

echo "Plotting RMSD figures..."
Rscript RMSDPlots.R $my_path $breakpoint_experiment $exp_num $chromosome