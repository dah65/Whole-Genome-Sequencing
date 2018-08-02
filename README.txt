#Next generation sequencing pipeline

#The purpose of this pipeline is to take raw sequence reads from a next generation sequencing platform 
#and perform trimming, mapping, read filtering, variant calling, and variant filtering. 

#Before you begin

#1) Upload all raw sequencing reads into the data/ directory. Reads from the same sample but different 
#sequencing lanes should be concatonated together before proceeding with the pipeline. 

#1) a. If your files are from multiple sequencing lanes and gunzipped (.gz extension), then you can use
#pbs_scripts/cat_reads.pbs to combine your files together for R1 and R2 separately. Make sure that your
#raw data is organized in the data/ directory accoriding to the .pbs script.

#2) Put your reference genome in the genomes/ directory and update the GENOME entry in the config.yaml file 

#2) a. If you only have a fasta file for your reference genome, you can use pbs_scripts/format_reference_genome.pbs
#but you will want to make sure that you replace p_humanus_reference.fa with the name of your reference genome.

#3) In the config.yaml file, update the SAMPLES entry to reflect all the sample names you would like included in the current analysis.

#4) Make sure the snakemake is installed in your environment. If not, run the following command:
#pip3 install --upgrade --ignore-installed --user snakemake

#5) To run snakemake, run pbs_scripts/snakemake_run.pbs, or pbs_scripts/snakemake_test.pbs to do a dry run.  
