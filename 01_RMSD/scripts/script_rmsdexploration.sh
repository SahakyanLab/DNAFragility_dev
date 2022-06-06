#!/bin/bash

my_path="$(pwd)/"

for kmer in 4 6 8
do
	Rscript RMSDExploration.R $my_path $kmer
done