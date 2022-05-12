#!/bin/bash

my_path="$(pwd)/"
action=$1
auto_fit=$2
upper_limit=$3

# Obtain heatmap tracts
# RMSD
for kmer in 4 6 8
do
    Rscript RMSDTracts.R $my_path $kmer $auto_fit
    Rscript KmertoneTractForML.R $my_path $kmer $action $upper_limit
done

# Kmertone
for kmer in 2 10 
do
    Rscript KmertoneTractForML.R $my_path $kmer $action $upper_limit
done