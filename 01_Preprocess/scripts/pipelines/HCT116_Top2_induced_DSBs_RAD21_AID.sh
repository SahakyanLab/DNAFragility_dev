#####################################################################################
# download sra data with 8 threads
parallel-fastq-dump --sra-id SRX5657096 --threads 8 --split-files --gzip

# to find the downloaded fastq files with sratoolkit
mv ~/ncbi/public/sra/* .

# rename downloaded fastq files
i=1
for file in *fastq*
do  
    mv $file reads_${i}.fastq.gz
    i=$((i+1))
done

i=1
for file in *fastq*
do  
    gunzip $file
    i=$((i+1))
done

# fastq trimming
mkdir fastq_trimmed
trim_galore -j 8 -q 20 --fastqc reads_1.fastq.gz -o fastq_trimmed -basename reads_1

# uses hg19 genome
# gunzip genome.fna.gz
path_to_ref="../../../data/ref/ref_index"

# build index
# $ bowtie-build reference_sequence.fasta index_name
# bowtie-build $ref genome.fa

# genome alignment
bowtie -l 50 -n 3 --best --strata --all -k 1 -t --sam ${path_to_ref}/genome.fa fastq_trimmed/reads_1_trimmed.fq.gz > output.sam.gz

# convert alignment output sam files to bam files and sort bam files 
samtools sort output.sam -o output_sorted.bam

# remove duplicates
samtools rmdup -s output_sorted.bam remove_duplicates.bam

# convert bam to bed files
bedtools bamtobed -i remove_duplicates.bam > remove_duplicates.bed

# sort bed file
sort -k 2,3 -u remove_duplicates.bed > output.bed

# # peak calling
# macs2 callpeak -t output_sorted.bam -f BAM -p 0.05 --keep-dup all -nomodel --outdir macs2_runs -nolambda

# fastq trimming
mkdir fastq_trimmed
trim_galore -j 8 -q 20 --fastqc reads_1.fastq.gz -o fastq_trimmed -basename reads_1

# genome alignment
bowtie -l 50 -n 3 --best --strata --all -k 1 -t --sam ${path_to_ref}/genome.fa fastq_trimmed/reads_1_trimmed.fq.gz > output.sam.gz

# convert alignment output sam files to bam files and sort bam files 
samtools sort output.sam.gz -o output_sorted.bam.gz

# remove duplicates
samtools markdup -r -s output_sorted.bam.gz remove_duplicates.bam.gz

# convert bam to bed files
bedtools bamtobed -i remove_duplicates.bam.gz > output.bed.gz