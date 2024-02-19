#!/bin/bash
mkdir -p ./bam/
mkdir -p ./filterbed/

for((i=1;i<=22;i++))
do
    wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr${i}-reali.bam.bai
    mv chr${i}-reali.bam.bai ./bam/chr${i}.bam.bai
    wget http://ftp.eva.mpg.de/neandertal/Chagyrskaya/BAM/chr${i}-reali.bam
    mv chr${i}-reali.bam ./bam/chr${i}.bam

    wget http://cdna.eva.mpg.de/neandertal/Chagyrskaya/FilterBed/chr${i}_mask.bed.gz
    mv chr${i}_mask.bed.gz ./filterbed/chr${i}_mask.bed.gz
    gunzip ./filterbed/chr${i}_mask.bed.gz      
done
