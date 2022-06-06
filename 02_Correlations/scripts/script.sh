#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path=$3
cores=$4
upper_limit=$5
first_idx=$6
breakpoint_path="../Raw_data/${breakpoint_experiment}/kmertone"

# create new folders to store calculations
mkdir -p ../data/{weight_factor,kmertone}/{$breakpoint_type,$breakpoint_experiment}
mkdir -p ../figures/{$breakpoint_type,$breakpoint_experiment}
mkdir -p ../data/seq_influence/

output_path="../02_Correlations/data/kmertone/${breakpoint_experiment}"

# Run kmertone software
# for kmer in 2 4 6 8
for kmer in 4
do	
	echo "Running kmertone enrichment/depletion analysis for kmer $kmer..."
	Rscript ../../Kmertone/Kmertone_run.R $my_path $breakpoint_path $kmer $ref_path $output_path $cores $first_idx
	
	echo "Obtaining weight factors for kmer $kmer..."
	Rscript ../lib/KmerWeightFactor.R $my_path $breakpoint_experiment $kmer
done

# # Generate plots
# for kmer in 2 4 6 8
# do
# 	Rscript KmertoneCorrelation.R $my_path $breakpoint_experiment $kmer $upper_limit

# 	for action in "ratio" "z-score"
# 	do
# 		Rscript KmertoneTractForML.R $my_path $breakpoint_experiment $kmer $upper_limit $action
# 	done
# done