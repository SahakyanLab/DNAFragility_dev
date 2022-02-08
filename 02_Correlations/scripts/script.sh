#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path=$3
cores=$4
exp_num=$5
upper_limit=$6
breakpoint_path="../../00_ReadsAlignment/data/kmertone/${breakpoint_experiment}_${exp_num}"
experiment=${breakpoint_experiment}_${exp_num}

# create new folders to store calculations
mkdir -p ../data/{weight_factor,kmertone}/{$breakpoint_type,$experiment}
mkdir -p ../figures/$breakpoint_type

output_path="../02_Correlations/data/kmertone/${breakpoint_experiment}_${exp_num}"

# Run kmertone software
for kmer in 2 4 6 8 10
do	
	echo "Running kmertone enrichment/depletion analysis for kmer $kmer..."
	Rscript ../../Kmertone/Kmertone_run.R $my_path $breakpoint_path $kmer $ref_path $output_path $cores
	
	echo "Obtaining weight factors for kmer $kmer..."
	Rscript ../lib/KmerWeightFactor.R $my_path $breakpoint_experiment $exp_num $kmer
done

# Generate plots
for kmer in 2 4 6 8
do
	Rscript KmertoneCorrelation.R $my_path $breakpoint_experiment $kmer $upper_limit

	for action in "ratio" "z-score"
	do
		Rscript KmertoneTractForML.R $my_path $breakpoint_experiment $kmer $upper_limit $action
	done
done