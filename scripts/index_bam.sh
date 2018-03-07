#!/bin/bash

# -------------------------------------------------------------------------------------- #
#Index new .bam file for aligned reads using samtools index
# -------------------------------------------------------------------------------------- #

#Check that the .bam file was passed as a parameter
USAGE="$0 aligned_reads.bam"
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

#Load required modules
module load samtools

#Move to correct directory
cd ./data

#Use samtools index to index your .bam with mapped reads
samtools index ./aligned_reads/$1 ./aligned_reads/${1}.bai

exit;
