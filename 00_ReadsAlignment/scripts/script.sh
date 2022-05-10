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

if [[ ! -f *.bam ]];
then
	if [ $nr_of_BAM_files -ge 1 ] && [ $nr_of_fasta_files -eq 0 ];
	then		
		echo "Extracting raw sequences from BAM file..."
		Rscript ../lib/ExtractBAMFiles.R $my_path $breakpoint_experiment
		BAM="TRUE"
		rm ../../Raw_data/$breakpoint_experiment/*.bam
	fi
fi

# align reads to reference sequence
for ((var=1; var<=$nr_of_fasta_files; var++))
do    
	echo "Counting number of reads..."
	fasta_lines="$(zgrep -Ec ">" ../../Raw_data/$breakpoint_experiment/chr$var.fasta.gz)"
	ind=$(((($fasta_lines+$interval-1)/$interval)-1))

	if [ $ind -lt 1 ]
	then
		ind=1
	fi

	if [ $nr_of_BAM_files -gt 0 ]
	then
		BAM="TRUE"
	else
		BAM="FALSE"
	fi

	echo "Aligning reads for chromosome $var..."
	Rscript ProcessingReads.R $my_path $breakpoint_experiment $var $ref_path $ind $fasta_lines $interval $BAM $alignment_strands
done

if [ $nr_of_fasta_files -eq 22 ]
then
	echo "Calculating levenshtein distances..."
	Rscript ../lib/CalcAverageLevenDist.R $my_path $breakpoint_experiment $cores
fi

# concatenate breakpoints into single file
echo "Concatenating breakpoints..."
Rscript ../lib/ConcatBreakpoints.R $my_path $breakpoint_experiment $cores $nr_of_fasta_files

# delete only when all csv files have been created
cd ../../Raw_data/$breakpoint_experiment/breakpoint_positions/
nr_of_files=$((ls *.csv | wc -l) 2> /dev/null)

if [[ ! -f *.csv ]]
then
	if [ $nr_of_files -eq 22 ]
	then
		rm -r */
	fi
fi

cd $my_path

if [ $nr_of_files -eq 22 ]
then
	echo "Obtaining breakpoint locations for kmertone..."
	Rscript ../lib/BreakpointsForKmertone.R $my_path $breakpoint_experiment $cores
fi