#!/bin/bash
mkdir -p ./bam/
mkdir -p ./filterbed/

wget http://cdna.eva.mpg.de/ust-ishim/BAM/Ust_Ishim.hg19_1000g.all.bam
mv Ust_Ishim.hg19_1000g.all.bam ./bam/Ust_Ishim.hg19_1000g.all.bam
samtools index ./bam/Ust_Ishim.hg19_1000g.all.bam
# parallel samtools index ::: ./bam/Ust_Ishim.hg19_1000g.all.bam

for((i=1;i<=22;i++))
do
    wget http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/Ust_Ishim/chr${i}_mask.bed.gz
    mv chr${i}_mask.bed.gz ./filterbed/chr${i}_mask.bed.gz
    gunzip ./filterbed/chr${i}_mask.bed.gz  
done
