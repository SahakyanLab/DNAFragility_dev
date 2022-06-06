#!/bin/bash

my_path="$(pwd)/"

# Process files
Rscript Extract_BSgenomes.R $pwd
Rscript Extract_savanna_elephant.R $pwd