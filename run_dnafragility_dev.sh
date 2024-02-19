#!/bin/bash

# preprocess files
cd ./01_Preprocess/scripts
bash submit.sh
cd ../../

# align reads
cd ./02_Alignment/scripts
bash submit.sh
cd ../../

# fit curves
cd ./03_FitCurves/scripts
bash submit.sh
cd ../../

# run kmertone
cd ./04_KmericAnalysis/scripts
bash submit.sh
cd ../../