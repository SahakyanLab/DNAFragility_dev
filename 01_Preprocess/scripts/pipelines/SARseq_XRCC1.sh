#####################################################################################
# download sra data with 8 threads
parallel-fastq-dump --sra-id SRX10151358 --threads 8 --split-files --gzip

# to find the downloaded fastq files with sratoolkit
mv ~/ncbi/public/sra/* .

# rename downloaded fastq files
i=1
for file in *fastq*
do  
    mv $file reads_${i}.fastq.gz
    i=$((i+1))
done

path_to_ref="../../../data/ref/ref_index"

# fastq trimming
mkdir fastq_trimmed
trim_galore -j 8 -q 20 --fastqc reads_1.fastq.gz -o fastq_trimmed -basename reads_1

# genome alignment
bowtie -l 50 -n 2 -m 1 -t --sam ${path_to_ref}/genome.fa fastq_trimmed/reads_1_trimmed.fq.gz > output.sam.gz

# convert alignment output sam files to bam files and sort bam files 
samtools sort output.sam.gz -o output_sorted.bam.gz

# remove duplicates
samtools markdup -r -s output_sorted.bam.gz remove_duplicates.bam.gz

# convert bam to bed files
bedtools bamtobed -i remove_duplicates.bam.gz > output.bed.gz