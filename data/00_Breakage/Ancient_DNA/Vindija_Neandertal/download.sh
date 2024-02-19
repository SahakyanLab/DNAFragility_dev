#!/bin/bash
mkdir -p ./bam/
mkdir -p ./filterbed/

for((i=1;i<=22;i++))
do
    wget http://ftp.eva.mpg.de/neandertal/Vindija/bam/Pruefer_etal_2017/Vindija33.19/Vi33.19.chr${i}.indel_realn.bam.bai
    mv Vi33.19.chr${i}.indel_realn.bam.bai ./bam/chr${i}_temp.bam.bai
    wget http://ftp.eva.mpg.de/neandertal/Vindija/bam/Pruefer_etal_2017/Vindija33.19/Vi33.19.chr${i}.indel_realn.bam
    mv Vi33.19.chr${i}.indel_realn.bam ./bam/chr${i}_temp.bam
    wget http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/Vindija33.19/chr${i}_mask.bed.gz
    mv chr${i}_mask.bed.gz ./filterbed/chr${i}_mask.bed.gz
    gunzip ./filterbed/chr${i}_mask.bed.gz
done