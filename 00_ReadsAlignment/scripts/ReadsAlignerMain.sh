#!/bin/bash

my_path="Documents/DPhil/SonicBreaks_dev/00_ReadsAlignment/scripts/"
cores=4
interval=5000000

for exp_num in 1
do
	reads_path="Simons_exp_${exp_num}"
	ref_path="Simons_exp"

	# create directories to save files
	mkdir ../data/reads/$reads_path/breakpoint_positions
	mkdir ../data/average_levdist/$reads_path
	mkdir ../../data/kmertone/$reads_path/breakpoints

	for var in {1..22}
	do
	    mkdir ../data/reads/$reads_path/breakpoint_positions/chr$var
	done

	# align reads to reference sequence
	for var in {1..22}
	do
	  echo "Starting read alignment for chromosome $var..."
	  
		fasta_lines="$(zgrep -Ec ">" ../data/reads/$reads_path/chr$var.fasta.gz)"
		ind=$(((($fasta_lines+$interval-1)/$interval)-1))
		
		for (( i=0; i<=$ind; i++ )) 
		do
			echo "Running read alignment script iteration $i of $ind..."
			Rscript ../lib/ReadsAligner.R $var $i $ref_path $reads_path $my_path $cores $fasta_lines $interval
		done
		
	done
	
	# obtain average levenshtein distance for all sequences
	mkdir ../figures/$reads_path/
	echo "Calculating levenshtein distances..."
	Rscript ../lib/LevenDistCalc.R $reads_path $my_path $cores
	echo "Calculated levenshtein distances!"

	# obtain files of breakpoint locations and chromosome number for kmertone
	echo "Obtaining breakpoint locations for kmertone..."
	Rscript ../lib/Kmertone_BreakPointLocations.R $reads_path $my_path $cores
	echo "Breakpoint locations for kmertone obtained!"
done
