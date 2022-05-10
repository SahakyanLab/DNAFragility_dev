#!/bin/bash

my_path="$(pwd)/"
cores=2
breakpoint_type="00-Ultrasonication"
exp="Simons_exp"
breakpoint_experiment="${breakpoint_type}/${exp}"
ref_path="Simons_exp"
chromosome=1
lower_limit=1
upper_limit=25

# create new folders to store calculations
mkdir -p ../data/{weight_factor,kmertone,overlap}/{$experiment,$breakpoint_type}

# obtain correlation matrix
for exp in $(seq $lower_limit $upper_limit)
do
  echo "Calculating overlap of exact breakage sites for experiment $exp..."
  Rscript ../lib/BreakpointOverlapCalc.R $my_path $breakpoint_experiment $ref_path $exp $lower_limit $upper_limit $cores
done

# obtain correlation plot from above matix
echo "Obtaining correlation matrix plot..."
Rscript ../lib/OverlapMatrix.R $my_path $experiment_folder $upper_limit