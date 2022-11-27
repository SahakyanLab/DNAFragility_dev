#####################################################################################
# download sra data with 8 threads
parallel-fastq-dump --sra-id SRR11680395 --threads 8 --split-files --gzip

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
trim_galore -j 8 -q 20 --fastqc reads_1.fastq -o fastq_trimmed -basename reads_1

# uses hg19 genome
gunzip genome.fna.gz
ref="genome.fna"

# build index
# $ bowtie-build reference_sequence.fasta index_name
bowtie2-build $ref genome.fa

# genome alignment
bowtie2 -N 0 -k 1 -q --local --fast-local -x genome.fa -U fastq_trimmed/reads_1.fq -S output.sam

# convert alignment output sam files to bam files and sort bam files 
samtools sort output.sam -o output_sorted.bam

# mark and remove duplicates
java -jar picard.jar MarkDuplicates I=output_sorted.bam O=output_sorted.marked.bam M=output_sorted.markdup.txt AS=true

# convert bam to bed files
bedtools bamtobed -i output_sorted.marked.bam > output_removed_duplicates.bed

# # peak calling
# macs2 callpeak -t output_sorted.bam -f BAM -p 0.05 --keep-dup all -nomodel --outdir macs2_runs -nolambda

path_to_ref="../../../data/ref/ref_index_bowtie2"

# fastq trimming
mkdir fastq_trimmed
trim_galore -j 8 -q 20 --fastqc reads_1.fastq.gz -o fastq_trimmed -basename reads_1

# genome alignment
bowtie2 -N 0 -k 1 -q --local --fast-local -x ${path_to_ref}/genome.fa -U reads_1_trimmed.fq.gz -S output.sam.gz

# convert alignment output sam files to bam files and sort bam files 
samtools sort output.sam.gz -o output_sorted.bam.gz

# remove duplicates
samtools markdup -r -s output_sorted.bam.gz remove_duplicates.bam.gz

# convert bam to bed files
bedtools bamtobed -i remove_duplicates.bam.gz > output.bed.gz