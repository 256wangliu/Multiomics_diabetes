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


module load bowtie/1.2.2
module load star/2.7.3a

mkdir -p $projPath2/fastq
cd $projPath2/fastq

sample_id='B13_D2 B13_D3 B14_D3 B14_D6 B15_D1 B15_D2 B16_D2 B16_D6 B6_D13 B6_D9 B7_D4 B7_D8 B8_D4 B8_D6 B13_R2 B13_R3 B14_R3 B14_R6 B15_R1 B15_R2 B16_R2 B16_R6 B6_R13 B6_R9 B7_R4 B7_R8 B8_R4 B8_R6'

for var in $sample_id

do
    histName=$var
cp ${projPath1}/rawdata/X202SC22082956-Z01-F001/usftp21.novogene.com/01.RawData/${histName}/${histName}_* ${projPath2}/fastq
cp ${projPath1}/rawdata/X202SC22082956-Z01-F002/usftp21.novogene.com/01.RawData/${histName}/${histName}_* ${projPath2}/fastq

mv ${projPath2}/fastq/${histName}_*_1.fq.gz ${projPath2}/fastq/${histName}_R1.fq.gz
mv ${projPath2}/fastq/${histName}_*_2.fq.gz ${projPath2}/fastq/${histName}_R2.fq.gz

done

for var in $sample_id
do
    histName=$var

# ## 1. fastqc
echo "1. FastQC"

mkdir -p ${projPath2}/fastqFileQC/${histName}
module load  fastqc/0.11.8

fastqc -o ${projPath2}/fastqFileQC/${histName} -f fastq ${projPath2}/fastq/${histName}_R1.fq.gz
fastqc -o ${projPath2}/fastqFileQC/${histName} -f fastq ${projPath2}/fastq/${histName}_R2.fq.gz

done 


for s in ${sample_id}
do

mkdir -p ${projPath2}/sam

p="$projPath/index/cell_id/cell_id" 

$projPath/Paired-Tag/reachtools/reachtools combine2 ${s}

zcat ${s}_combined.fq.gz | bowtie ${p} - --norc -m 1 -v 1 -S ${projPath2}/sam/${s}_BC.sam

#### This step convert to Celluar Barcode mapped reads to fastq files.
$projPath/Paired-Tag/reachtools/reachtools convert2 ${projPath2}/sam/${s}_BC.sam  
done


