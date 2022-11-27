#!/bin/bash

my_path="$(pwd)/"
breakpoint_type=$1
exp=$2
file_format=$3
start_idx=$4
end_idx=$5
alignment_strands=$6
fastq_processed=$7

# cthreepo -i olambda_summits.bed -if rs -it uc -f bed -m h37 -ku -o output.bed
# cthreepo -i output.bed -if rs -it uc -f bed -m h37 -ku -o out.bed

if [ "${fastq_processed}" == "TRUE" ]
then
    path_to_file="../../Raw_data/${breakpoint_type}/${exp}"
    if [ -f ${path_to_file}/output.bed.gz ]
    then
        mv ${path_to_file}/output.bed.gz ${path_to_file}/output.bed
        check=$(head -1 output.bed | awk '{print $1}')
        if [ ${check} != "chr1" ]
        then
            echo "Converting coordinates..."
            cthreepo -i ${path_to_file}/output.bed -if rs -it uc -f bed -m h37 -ku -o ${path_to_file}/out.bed
            mv ${path_to_file}/out.bed ../../Raw_data/${breakpoint_type}/${exp}/output.bed
        fi
    fi
fi

# Get genomes
Rscript Process.R $my_path "FALSE"

# Run all pipelines
cd ./pipelines
for file in ./*.sh
do
    bash $file
done
cd ../

# Process files
Rscript Process.R $my_path "TRUE"