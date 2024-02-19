#!/usr/bin/bash

for folder in ./*
do
    cd "$folder"
    echo "Processing files in $folder"
    pwd="$(pwd)/"
    Rscript Process.R $pwd
    cd ../
done