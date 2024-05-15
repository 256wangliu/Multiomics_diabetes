#!/bin/bash
#$ -pe threaded 8
#$ -N Paired-Tag
#$ -q 4-day
#$ -j y
#$ -o /research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/tmp
#$ -e /research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/tmp
#$ -M Wang.Liu@mayo.edu
#$ -m abe
#$ -l h_vmem=8G
#$ -notify

TMPDIR= "/research/labs/bme/weiz/m238739/Paired-Tag/"
projPath="/research/labs/bme/weiz/m238739/Paired-Tag/"
projPath1="/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/"
projPath2="/research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse"

#####################
######RNA
mm10="$projPath/index/star_GRCm38"
mm10_rna="$projPath/index/star_GRCm38/mm10.big.txt"


mkdir -p $projPath2/bam/RNA
sample_id='B13_R2 B13_R3 B14_R3 B14_R6 B15_R1 B15_R2 B16_R2 B16_R6 B6_R13 B6_R9 B7_R4 B7_R8 B8_R4 B8_R6'

module load cutadapt/2.8

module load samtools/1.10

module load star/2.7.3a

for s in ${sample_id}

do

cd $projPath2/fastq

$projPath/TrimGalore-0.6.6/trim_galore $projPath2/sam/${s}_BC_cov.fq.gz

$projPath/TrimGalore-0.6.6/trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${s}_BC_cov_trimmed.fq.gz ### trim oligo-dT primer

$projPath/TrimGalore-0.6.6/trim_galore -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${s}_BC_cov_trimmed_trimmed.fq.gz ## trim N6 primer

STAR  --runThreadN 6 --genomeDir ${mm10} --readFilesIn ${s}_BC_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix $projPath2/bam/RNA/${s}_mm10_ --outSAMtype BAM Unsorted

samtools view -h -F 256 $projPath2/bam/RNA/${s}_mm10_Aligned.out.bam -b > $projPath2/bam/RNA/${s}\_clean.bam

samtools view -h $projPath2/bam/RNA/${s}\_clean.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $projPath2/bam/RNA/${s}\_clean_tmp.bam

samtools sort $projPath2/bam/RNA/${s}\_clean_tmp.bam -o $projPath2/bam/RNA/${s}_sorted.bam

cd $projPath2/bam/RNA

$projPath/Paired-Tag/reachtools/reachtools rmdup2 $projPath2/bam/RNA/${s}_sorted.bam

$projPath/Paired-Tag/reachtools/reachtools bam2Mtx2 $projPath2/bam/RNA/${s}_sorted_rmdup.bam ${mm10_rna}

done
