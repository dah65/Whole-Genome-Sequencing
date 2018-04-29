#!/bin/bash

# -------------------------------------------------------------------------------------- #
#Trim raw paried end sequencing reads using trimmomatic
# -------------------------------------------------------------------------------------- #

#Check that forward and reverse reads have been passed as parameters
USAGE="$0 <Sample ID>"
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1; 
fi

#Load trimmomatic
module load trimmomatic

#Move into directory with raw data
cd ./data

#Define parameters as forward and reverse reads
READ1=${1}_R1.fastq
READ2=${1}_R2.fastq

#Set path to adapter sequences
ADAPTERS=/storage/work/cxb585/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa

#Run trimmomatic
java -jar $TRIMMOMATIC PE -phred33 $READ1 $READ2 ./trimmed_reads/${READ1}_paired.fastq.gz ./trimmed_reads/${READ1}_unpaired.fastq.gz ./trimmed_reads/${READ2}_paired.fastq.gz ./trimmed_reads/${READ2}_unpaired.fastq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

exit;
