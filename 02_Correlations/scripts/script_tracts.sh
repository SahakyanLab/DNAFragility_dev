#!/bin/bash

my_path="$(pwd)/"
action=$1
auto_fit=$2
upper_limit=$3
iter=$4

# Obtain heatmap tracts
# RMSD
for kmer in 4 6 8
do
    Rscript RMSDTracts.R $my_path $kmer $auto_fit $iter
done

# Kmertone
if [ "${auto_fit}" == "FALSE" ]
then
    for kmer in 2 4 6 8
    do
        Rscript KmertoneTractForML.R $my_path $kmer $action $upper_limit
    done
fi