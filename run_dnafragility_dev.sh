#!/bin/bash

# preprocess files
cd ./01_Process
bash submit.sh
cd ../

# align reads
cd ./02_Alignment
bash submit.sh
cd ../

# fit curves
cd ./03_FitCurves
bash submit.sh
cd ../

# run kmertone
cd ./04_KmerAnalysis
bash submit.sh
cd ../