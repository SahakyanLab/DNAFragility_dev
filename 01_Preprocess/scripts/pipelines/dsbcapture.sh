#####################################################################################
# download sra data with 8 threads
# NHEK_break_seq_rep1
# Rename to: BREAK_rep1_R1.fastq and BREAK_rep1_R2.fastq
# parallel-fastq-dump --sra-id SRR3182242 --threads 8 --split-files --gzip 
# fasterq-dump --split-files SRR3182242

# NHEK_break_seq_rep2
# Rename to: BREAK_rep2_R1.fastq and BREAK_rep2_R2.fastq
# parallel-fastq-dump --sra-id SRR3182246 --threads 8 --split-files --gzip
# fasterq-dump --split-files SRR3182246

# NHEK_bless_seq_rep1
# Rename to: BLESS_rep1_R1.fastq and BLESS_rep1_R2.fastq
# parallel-fastq-dump --sra-id SRR3182247 --threads 8 --split-files --gzip 
# fasterq-dump --split-files SRR3182247

# NHEK_bless_seq_rep2 
# Rename to: BLESS_rep2_R1.fastq and BLESS_rep2_R2.fastq
# parallel-fastq-dump --sra-id SRR3182252 --threads 8 --split-files --gzip 
# fasterq-dump --split-files SRR3182252

# AID-DlvA_AsiSI_restriction_enzyme
# Rename to: AID-DlvA_AsiSI_R1.fastq and AID-DlvA_AsiSI_R2.fastq
# fasterq-dump --split-files SRR3182238

# EcoRV_restriction_enzyme
# Rename to: EcoRV_restriction_enzyme.fastq
fasterq-dump --split-files SRR3182237

#####################################################################################
#####################################################################################
######################################## TEST #######################################
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

# conda create -n samtools samtools
# conda activate samtools

# ref="/project/sahakyanlab/ppflugha/03_Breakpoints/Raw_data/ref/ref_index/genome.fa"
ref="../../ref/ref_index_bwa/genome.fa"
whitelist="wgEncodeDukeMapabilityRegionsExcludable.bed"

# BREAK-seq (DSBCapture)
echo "Trim adapters and filter bases based on quality"
for fq1 in BREAK*R1*fastq; 
do
  fq2=${fq1/_R1/_R2}
  out=`basename $fq1 _R1.fastq`
  out1=${out}.R1.fq
  out2=${out}.R2.fq
  cutadapt -j 2 -O 3 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $out1 -p $out2 $fq1 $fq2
done

echo "Align with BWA MEM"
for f in *R1.fq;
do
    bwa mem -M -t 10 $ref $f ${f/R1.fq/R2.fq} | samtools view -@ 10 -L $whitelist -F 2820 -Sb - > ${f%%.R1.fq}.bam
done

echo "Sorting and merging"
for f in *.bam;
do 
    samtools sort -@ 10 $f -o ${f%%.bam}.sort
done
  
samtools merge -@ 10 BREAK-primary-AD12_S1.sort.bam BREAK*sort

# brew install picard-tools
# picard MarkDuplicates
echo "Mark duplicate reads"
for f in *sort.bam;
do
  picard MarkDuplicates I="$f" O="${f%%.bam}.marked.bam" M="${f%%.bam}.markdup.txt" AS=true
done

# Perform quality filter (q>10), n-sorting, 
# Fixmates to remove unpaired reads
# Get only first read (first read extraction only for DSBCapture)
echo "Perform quality filter"
for f in *marked.bam;
do
  f2=${f/sort/nsort}
  f3=${f2%%.bam}
  samtools sort -@ 10 -n $f -o ${f2%%.bam};
  samtools view -b -@ 10 -F 1024 -q 10 ${f2%%.bam} | samtools fixmate - ${f3}.cleaned.fixed.bam;
  samtools view -f 65 -b ${f3}.cleaned.fixed.bam > ${f3}.cleaned.fixed.first.bam;
  samtools sort -@ 10 ${f3}.cleaned.fixed.first.bam -o ${f3}.cleaned.fixed.first.sort.bam;
  samtools index ${f3}.cleaned.fixed.first.sort.bam
done
######################################## TEST #######################################
#####################################################################################
#####################################################################################

# ref="genome.fa"
# echo "Creating BWA-MEM index..."
# bwa index $ref

# ref="/project/sahakyanlab/ppflugha/03_Breakpoints/Raw_data/ref/ref_index/genome.fa"
# whitelist="wgEncodeDukeMapabilityRegionsExcludable.bed"

# # BREAK-seq (DSBCapture)
# echo "trim adapters and perform alignment discarding unaligned..."
# for fq1 in BREAK*R1*fastq
# do
#   fq2=${fq1/R1/R2}
#   out=`basename $fq1 _R1.fastq`
#   out1=${out}.R1.fq
#   out2=${out}.R2.fq
#   cutadapt -O 3 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o $out1 -p $out2 $fq1 $fq2
# done

# echo "Now bwa..."
# for f in *R1.fq;
# do
#     bwa mem -M -t 10 $ref $f ${f/R1.fq/R2.fq} | samtools view -@ 10 -L $whitelist -F2820 -Sb - > ${f%%.R1.fq}.bam
# done

# # BLESS
# echo "trim adapters and keep only read with both linkers..."
# for fq1 in BLESS*R1*; 
# do
#   fq2=${fq1/_R1_/_R2_}
#   adir="./${fq1%%.fastq}" 
#   mkdir $adir
#   mkdir fastq_trimmed
#   cutadapt -m 10 -O 5 -g TCGAGGTAGTA -G TCGAGACGACG -o ${adir}/tmp_out1.fq -p ${adir}/tmp_out2.fq $fq1 $fq2 --discard-untrimmed;
#   cutadapt -m 10 -O 5 -G TCGAGGTAGTA -g TCGAGACGACG -o ${adir}/tmp_out3.fq -p ${adir}/tmp_out4.fq $fq1 $fq2 --discard-untrimmed;
#   cat ${adir}/tmp_out1.fq ${adir}/tmp_out3.fq > ${adir}/tmp_out13.fq;
#   cat ${adir}/tmp_out2.fq ${adir}/tmp_out4.fq > ${adir}/tmp_out24.fq;
#   cutadapt -m 10 -O 5 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o fastq_trimmed/${fq1%%.fastq}.trimmed.fq -p fastq_trimmed/${fq2%%.fastq}.trimmed.fq ${adir}/tmp_out13.fq ${adir}/tmp_out24.fq
# done

# for f in `find . -name 'tmp_out13.fq'`
# do
#   f1=${f##./}
#   f2=${f1////_}
#   mv $f $f2  
# done

# for f in `find . -name 'tmp_out24.fq'`
# do
#   f1=${f##./}
#   f2=${f1////_}
#   mv $f $f2  
# done

# echo "Now bwa..."
# for f in *tmp_out13.fq;
# do
#   bwa mem -M -t 10 $ref $f ${f/out13/out24} | samtools view -@ 10 -L $whitelist -F2820 -Sb - > ${f%%_tmp*}.bam
# done

# echo "For both BLESS and DSBCapture post-alignment processing..."
# for f in *.bam; do echo ${f%%.bam}.sort ; done
# do 
#     samtools sort -@ 10 $f -o ${f%%.bam}.sort
# done
  
# samtools merge -@ 10 BREAK-primary-AD12_S1.sort.bam BREAK*sort
# samtools merge -@ 10 BLESS-primary-AD16_S4_R1_001.sort.bam BLESS*sort

# echo "Perform markduplicates..."
# for f in *sort.bam; 
# do
#   MarkDuplicates I="$f" O="${f%%.bam}.marked.bam" M="${f%%.bam}.markdup.txt" AS=true
# done

# ############################################################################################################################
# # (q>10), n-sorting, fixmates to remove unpaired reads and get only first read (first read extraction only for DSBCapture) #
# ############################################################################################################################ 
# echo "perform quality filter..."
# for f in *marked.bam;
# do
#   f2=${f/sort/nsort}
#   f3=${f2%%.bam}
#   samtools sort -@ 10 -n $f -o ${f2%%.bam};
#   samtools view -b -@ 10 -F 1024 -q 10 ${f2%%.bam} | samtools fixmate - ${f3}.cleaned.fixed.bam;
#   samtools view -f 65 -b ${f3}.cleaned.fixed.bam > ${f3}.cleaned.fixed.first.bam;
#   samtools sort -@ 10 ${f3}.cleaned.fixed.first.bam -o ${f3}.cleaned.fixed.first.sort;
#   samtools index ${f3}.cleaned.fixed.first.sort
# done

# echo "instead of extracting first read, extract proximal linker for BLESS..."
# for f in  BLESS*R1*fastq
# do
#   cutadapt -m 10 -O 5 -g TCGAGGTAGTA -o ${f%%.fastq}.id1.fq $f --discard-untrimmed
# done

# for f in  BLESS*R2*fastq
# do
#   cutadapt -m 10 -O 5 -g TCGAGGTAGTA -o ${f%%.fastq}.id2.fq $f --discard-untrimmed
# done

# for f in *R1*id1*
# do
#   grep -i '^@NS' $f > ${f%%.fq}.txt
# done

# for f in *R2*id2*
# do
#   grep -i '^@NS' $f > ${f%%.fq}.txt
# done

# mkdir split
# for f in *id1*txt
# do
#   split -l 1000000 -d $f ${f%%.txt}_split_
# done

# for f in *id2*txt
# do
#   split -l 1000000 -d $f ${f%%.txt}_split_
# done

# # cd /lustre/sblab/marsic01/repository/NextSeq/BREAK_BLESS_primary/n_2/split/rep_1
# # samtools view -h -f 65 BLESS_primary.SLX-10843.nsort.marked.cleaned.fixed.bam > BLESS_primary.SLX-10843.nsort.marked.cleaned.fixed.R1.sam
# # samtools view -h -f 129 BLESS_primary.SLX-10843.nsort.marked.cleaned.fixed.bam > BLESS_primary.SLX-10843.nsort.marked.cleaned.fixed.R2.sam

# # for f in *split*
# # do
# #    bsub "cut -f 1 -d ' '  $f | sed 's/@//' > ${f}.txt"
# # done

# # for f in *id1*txt
# # do
# #   bsub "split -l 1000000 -d $f ${f%%.txt}_split_"
# # done

# # cd /lustre/sblab/marsic01/repository/NextSeq/BREAK_BLESS_primary/n_1/split
# # for f in *id1_*txt
# # do
# #   bsub -R "rusage[mem=8000]" "LC_ALL=C fgrep -f $f BLESS-primary-AD16_S4_R1_001.nsort.marked.cleaned.fixed.R1.sam > ${f%%.txt}.R1.proximal.sam" 
# # done

# # for f in *id2_*txt
# # do
# #   bsub -R "rusage[mem=8000]" "LC_ALL=C fgrep -f $f BLESS-primary-AD16_S4_R1_001.nsort.marked.cleaned.fixed.R2.sam > ${f%%.txt}.R2.proximal.sam" 
# # done

# # bsub "head -n 28 BLESS-primary-AD16_S4_R1_001.nsort.marked.cleaned.fixed.R1.sam | cat - *R1.proximal.sam > BLESS-primary-AD16_S4_R1_001.nsort.marked.cleaned.fixed.R1.proximal.sam
# # bsub "head -n 28 BLESS-primary-AD16_S4_R1_001.nsort.marked.cleaned.fixed.R2.sam | cat - *R2.proximal.sam > BLESS-primary-AD16_S4_R1_001.nsort.marked.cleaned.fixed.R2.proximal.sam"
 
# # cd /lustre/sblab/marsic01/repository/NextSeq/BREAK_BLESS_primary/n_2/split/rep_1
# # for f in *id1_*txt
# # do
# #   bsub -R "rusage[mem=8000]" "LC_ALL=C fgrep -f $f BLESS_primary.SLX-10843.nsort.marked.cleaned.fixed.R1.sam > ${f%%.txt}.R1.proximal.sam" 
# # done
# # for f in *id2_*txt
# # do
# #   bsub -R "rusage[mem=8000]" "LC_ALL=C fgrep -f $f BLESS_primary.SLX-10843.nsort.marked.cleaned.fixed.R2.sam > ${f%%.txt}.R2.proximal.sam" 
# # done

# # # =======================
# # # peak calling
# # # =======================

# for f in *first.sort
# do
#   macs2 callpeak --keep-dup all -t $f -n ${f}.first.sort.NOdups.q005
# done