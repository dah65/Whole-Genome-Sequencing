#!/bin/bash

# -------------------------------------------------------------------------------------- #
#Convert aligned .sam file to .bam
# -------------------------------------------------------------------------------------- #

#Check that reference genome and .sam file have been passed as parameters
USAGE="$0 reference_genome.fa aligned.sam"
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

#Load required modules
module load samtools

#Move to working directory
cd ./data

#Set parameters
GENOME_INDEX=../genomes/${1}.fai
ALIGNED_READS=./aligned_reads/$2

#Run samtools view to convert .sam to .bam
samtools view -S -b -t $GENOME_INDEX -o ${ALIGNED_READS}.bam $ALIGNED_READS

exit;
