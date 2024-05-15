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

###############################
module load cutadapt/2.8
module load bowtie2/2.3.3.1


sample_id='B13_D2 B13_D3 B14_D3 B14_D6 B15_D1 B15_D2 B16_D2 B16_D6 B6_D13 B6_D9 B7_D4 B7_D8 B8_D4 B8_D6'

module load samtools/1.10

mkdir -p $projPath2/sam
mkdir -p $projPath2/bam


for s in ${sample_id}
do

module load python/3.9.2
mm10="$projPath/index/bowtie2_GRCm38/mouse"
mm10_5k="$projPath/index/bowtie2_GRCm38/mm10.bin5k.txt"

cd $projPath2/fastq

$projPath/TrimGalore-0.6.6/trim_galore $projPath2/sam/${s}_BC_cov.fq.gz

bowtie2 -x ${mm10} -U ${s}_BC_cov_trimmed.fq.gz --no-unal -p 8 -S $projPath2/sam/${s}_mm10.sam

samtools view -bS $projPath2/sam/${s}_mm10.sam > $projPath2/bam/${s}_mm10.bam

samtools view -h $projPath2/bam/${s}_mm10.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $projPath2/bam/${s}_mm10_tmp.bam

samtools sort -o  $projPath2/bam/${s}_mm10_sorted.bam  $projPath2/bam/${s}_mm10_tmp.bam

samtools index $projPath2/bam/${s}_mm10_sorted.bam

cd $projPath2/bam

$projPath/Paired-Tag/reachtools/reachtools rmdup2 $projPath2/bam/${s}_mm10_sorted.bam

samtools index $projPath2/bam/${s}_mm10_sorted_rmdup.bam

python $projPath/Paired-Tag/remove_pileup/count_pileups.py $projPath2/bam/${s}_mm10_sorted_rmdup.bam $projPath2/bam/${s}.filtered.pos_strand.tally.txt

python $projPath/Paired-Tag/remove_pileup/remove_pileups.py $projPath2/bam/${s}_mm10_sorted_rmdup.bam $projPath2/bam/${s}.filtered.pos_strand.tally.txt $projPath2/bam/${s}_mm10_filtered_10.bam 10

$projPath/Paired-Tag/reachtools/reachtools bam2Mtx2 $projPath2/bam/${s}_mm10_filtered_10.bam ${mm10_5k}

done


