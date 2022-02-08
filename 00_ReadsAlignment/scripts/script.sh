#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path=$3
cores=$4
interval=$5
alignment_strands=$6

for exp_num in 1
do
	reads_path="${breakpoint_experiment}_${exp_num}"

	# create directories to save files
	mkdir -p ../data/{reads,average_levdist}/$breakpoint_type
	mkdir -p ../data/kmertone/$breakpoint_type
	mkdir -p ../data/reads/$reads_path/breakpoint_positions/chr{1..22}
	mkdir -p ../data/reads/$reads_path/breakpoint_positions/concat
	mkdir -p ../figures/{$breakpoint_type,$reads_path}

	# check if any BAM files present
	cd ../data/reads/$reads_path/
	nr_of_files=$((ls *.bam | wc -l) 2> /dev/null)

	if [[ ! -f *.bam ]];
	then
		if [ $nr_of_files -gt 1 ];
		then		
			echo "Extracting raw sequences from BAM file..."
			Rscript ../lib/ExtractBAMFiles.R $my_path $breakpoint_experiment $exp_num
		fi
	fi

	cd ../../../../scripts
	
	# align reads to reference sequence
	for var in {1..22}
	do    
		echo "Counting number of reads..."
		fasta_lines="$(zgrep -Ec ">" ../data/reads/$reads_path/chr$var.fasta.gz)"
		ind=$(((($fasta_lines+$interval-1)/$interval)-1))

		echo "Aligning reads for chromosome $var..."
		Rscript ProcessingReads.R $my_path $breakpoint_experiment $exp_num $var $ref_path $ind $fasta_lines $interval $BAM $alignment_strands
	done
	
	echo "Calculating levenshtein distances..."
	Rscript ../lib/CalcAverageLevenDist.R $my_path $breakpoint_experiment $exp_num $cores

	concatenate breakpoints into single file
	echo "Concatenating breakpoints..."
	Rscript ../lib/ConcatBreakpoints.R $my_path $breakpoint_experiment $exp_num $cores
	
	# delete only when all csv files have been created
	cd ../data/reads/$reads_path/breakpoint_positions/
	nr_of_files=$((ls *.csv | wc -l) 2> /dev/null)

	if [[ ! -f *.csv ]];
	then
		if [ $nr_of_files -eq 22 ];
		then
			rm -r */
		fi
	fi

	cd ../../../../../scripts

	echo "Obtaining breakpoint locations for kmertone..."
	Rscript ../lib/BreakpointsForKmertone.R $my_path $breakpoint_experiment $exp_num $cores
done