#!/bin/bash
mkdir -p ./bam/
mkdir -p ./filterbed/

wget http://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam
mv T_hg19_1000g.bam ./bam/T_hg19_1000g.bam
wget http://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam.bai
mv T_hg19_1000g.bam.bai ./bam/T_hg19_1000g.bam.bai

for((i=1;i<=22;i++))
do
    wget http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/Denisova/chr${i}_mask.bed.gz
    mv chr${i}_mask.bed.gz ./filterbed/chr${i}_mask.bed.gz
    gunzip ./filterbed/chr${i}_mask.bed.gz  
done