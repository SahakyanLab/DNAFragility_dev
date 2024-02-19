#!/bin/bash

pwd="$(pwd)/"

# Get genomes
Rscript Process.R $pwd "FALSE"

# Process files
Rscript Process.R $pwd "TRUE"