#!/bin/bash

my_path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/scripts/ReadsAlignment/"

for exp_num in 6
do
	pre_reads_path="Simons_exp_$((exp_num-1))"
	reads_path="Simons_exp_${exp_num}"
	ref_path="Simons_exp"
	two_mer_ref="Simons_exp_ref"

	# unzip raw read files
	echo "Unzipping raw read files..."
	gunzip ../../data/reads/$reads_path/*
	echo "Raw read files unzipped!"

	# create directories to save files
	mkdir ../../data/reads/$reads_path/breakpoint_positions
	mkdir ../../data/reads/$reads_path/fasta_chunks
	mkdir ../../data/reads/two_mers/$reads_path

	for var in {1..22}
	do
	    mkdir ../../data/reads/$reads_path/breakpoint_positions/chr$var
	done

	# align reads to reference sequence 
	for var in {1..22}
	do
		echo "Splitting fasta file into chunks of 1M lines..."
		awk -v size=1000000 -v pre=chunks -v pad=1 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' ../../data/reads/$reads_path/chr$var.fasta
		echo "Fasta file split!"
		mv chunks* ../../data/reads/$reads_path/fasta_chunks/

		for ind in 1 2
		do
			echo "Running python script..."
			python3 ReadsAligner.py $var $ind $ref_path $reads_path $my_path
		done

		# remove split fasta files
		rm ../../data/reads/$reads_path/fasta_chunks/*
		echo "Job done with chr$var"
	done

	# obtain average levenshtein distance for all sequences
	mkdir ../../figures/$reads_path/
	echo "Calculating levenshtein distances..."
	Rscript ../../lib/LevenDistCalc.R $reads_path $my_path
	echo "Calculated levenshtein distances!"

	# obtain normalised breakpoint data
	echo "Calculating normalised breakpoint frequencies..."
	Rscript ../../lib/BreakPointFrequency.R $reads_path $ref_path $my_path	
	echo "Calculated normalised breakpoint frequencies!"

	# obtain files of breakpoint locations and chromosome number for kmertone
	echo "Obtaining breakpoint locations for kmertone..."
	Rscript ../../lib/Kmertone_BreakPointLocations.R $reads_path $my_path
	echo "Breakpoint locations for kmertone obtained!"

	# zip read files to save hard disk space
	echo "Zipping raw reads files..."
	gzip ../../data/reads/$reads_path/*
	echo "Raw read files zipped!"
done
