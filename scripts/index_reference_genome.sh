#`/bin/bash

# -------------------------------------------------- #
#Index the reference genome using BWA and samtools
# -------------------------------------------------- #

#Check that genome was passed as a parameter
USAGE="$0 reference_genome.fa"
if [ -z $1 ]; then
	echo "ERROR $USAGE";
	exit 1;
fi

#Move to directory containing reference genome
cd ./genomes

#Load required modules
module load bwa samtools

#Index the genome using BWA
bwa index $1

#Index the genome using samtools
samtools faidx $1

exit;
