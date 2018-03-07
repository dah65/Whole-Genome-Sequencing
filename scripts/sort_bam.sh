#!/bin/bash

# -------------------------------------------------------------------------------------- #
#Sort new .bam file for aligned reads using samtools sort
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

#Use samtools sort to sort your .bam with mapped reads by location in genome
samtools sort ./aligned_reads/$1 -o ./aligned_reads/${1}.sorted.bam --threads 8

exit;
