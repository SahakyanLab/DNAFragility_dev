#!/bin/bash

pwd="$(pwd)/"

# Process files
for kmer in 6 8
do
    Rscript Process.R $pwd $kmer "zscore"
    Rscript Process.R $pwd $kmer "ratio"
done

# analysing results
Rscript 00_Extend_QueryTable.R $pwd
Rscript 01_Clustering_maps.R $pwd
Rscript 02_Breakage_Insights.R $pwd 8 "zscore" "FALSE"
Rscript 03_ExactBreaks.R $pwd
Rscript 03_ExactBreaks_and_KmerCorr_Plot.R $pwd