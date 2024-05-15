#!/bin/bash
#$ -pe threaded 8
#$ -N Paired-Tag
#$ -q 1-day
#$ -j y
#$ -o /research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/tmp
#$ -e /research/labs/bme/weiz/m238739/Paired-Tag/data/data16/mouse/tmp
#$ -M Wang.Liu@mayo.edu
#$ -m abe
#$ -l h_vmem=8G
#$ -notify

projPath="/research/labs/bme/weiz/m238739/Paired-Tag/"
sample_id='B13'
for s in ${sample_id}
do
cd $projPath/data/data16/mouse/bam

$projPath/Paired-Tag/perlscripts/merge_mtx.pl merge_list_DNA_${s}.txt

cd $projPath/data/data16/mouse/bam/RNA/

$projPath/Paired-Tag/perlscripts/merge_mtx.pl merge_list_RNA_${s}.txt

done












