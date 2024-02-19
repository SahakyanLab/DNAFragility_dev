#!/bin/bash
mkdir -p ./bam/
mkdir -p ./filterbed/

for((i=1;i<=22;i++))
do
    wget http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/AltaiNea.hg19_1000g.${i}.dq.bam.bai
    mv AltaiNea.hg19_1000g.${i}.dq.bam.bai chr${i}.bam.bai
    wget http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/bam/AltaiNea.hg19_1000g.${i}.dq.bam
    mv AltaiNea.hg19_1000g.${i}.dq.bam chr${i}.bam
    wget http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/Altai/chr${i}_mask.bed.gz
    mv chr${i}_mask.bed.gz ./filterbed/chr${i}_mask.bed.gz
    gunzip ./filterbed/chr${i}_mask.bed.gz  
done
