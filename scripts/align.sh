#!/bin/bash

# -------------------------------------------------------------------------------------- #
#Align reads to reference genome using BWA
# -------------------------------------------------------------------------------------- #

#Check that genome and reads have been passed as parameters
USAGE="$0 reference_genome.fa trimmed_R1.fa trimmed_R2.fa"
if [ -z $3 ]; then
	echo "ERROR $USAGE";
	exit 1;
fi

#Move to working directory
cd ./data

#Load required modules
module load bwa

#Name parameters
IND_ID=CAH9

#Run BWA mem
bwa mem -t 8 ../genomes/${1} trimmed_reads/${2} trimmed_reads/${3} > aligned_reads/${IND_ID}.PE.bwa.${1}.sam

exit;
