#!/bin/bash

pwd="$(pwd)/"

# Process files
for kmer in 6 8
do
    Rscript Process.R $pwd $kmer
done

# analysing results
Rscript 01_Clustering_curves.R $pwd
Rscript 04_RSS_per_curve.R $pwd

# if exists
if [ -d "../../../DNAFragility/COSMIC/data/annotations/" ]
then
    Rscript 05_Breaks_at_genic_sites.R $pwd
    Rscript 06_RMSD_without_promoters.R $pwd
fi

Rscript 07_1_RMSD_control.R $pwd
Rscript 08_RMSD_all_chr.R $pwd
Rscript 09_DNAShape_comparisons.R $pwd
Rscript 11_Calc_RMSD_Diff.R $pwd