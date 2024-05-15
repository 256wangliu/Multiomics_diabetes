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

cd $projPath/data/data16/mouse/bam
sample_id='B13_D2 B13_D3'

for s in ${sample_id}
do

metadata="$projPath/data/data16/mouse/bam/${s}_DNA_RNA_filt.xls"
pre_mtx="$projPath/data/data16/mouse/bam/${s}_mm10_filtered_10_mtx2"
prefix="${s}"

export metadata
export pre_mtx
export prefix
 
$projPath/Paired-Tag/perlscripts/filt_mtx.pl 
done


cd $projPath/data/data16/mouse/bam/RNA
sample_id='B13_D2:B13_R2 B13_D3:B13_R3'

for s in ${sample_id}
do

first=$(echo $s | cut -d ":" -f 1)
second=$(echo $s | cut -d ":" -f 2)

metadata="$projPath/data/data16/mouse/bam/${first}_DNA_RNA_filt.xls"
pre_mtx="$projPath/data/data16/mouse/bam/RNA/${second}_sorted_rmdup_mtx2"
prefix="${second}"

export metadata
export pre_mtx
export prefix

$projPath/Paired-Tag/perlscripts/filt_mtx.pl

done
