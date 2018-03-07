configfile: "./config.yaml"

rule load_modules:
  input:
    expand("{module}", module=config["MODULES"])
  shell:
    "module load {input}"

rule trim_reads:
  input:
    expand("data/{sample}_R1.fastq", sample=config["SAMPLES"]),
    expand("data/{sample}_R2.fastq", sample=config["SAMPLES"])
  output:
    "data/trimmed_reads/{sample}_R1_paired.fastq.gz",
    "data/trimmed_reads/{sample}_R1_unpaired.fastq.gz",
    "data/trimmed_reads/{sample}_R2_paired.fastq.gz",
    "data/trimmed_reads/{sample}_R2_unpaired.fastq.gz"
  log: "logs/trim/{sample}.trim.log"
  shell:
    "java -jar $TRIMMOMATIC PE -phred33 {input} {output} ILLUMINACLIP:{config[ADAPTERS]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    
rule bwa_map_and_samtools_make_bam:
  input:
    "genomes/p_humanus_reference.fa",
    expand("data/trimmed_reads/{sample}_R1_paired.fastq.gz", sample=config["SAMPLES"]),
    expand("data/trimmed_reads/{sample}_R2_paired.fastq.gz", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.bam"
  log: "logs/{sample}.map_and_bam.log"
  shell:
    "bwa mem -t 8 {input} | samtools view -Sb > {output}"

rule sort_bam:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.bam"
  log: "logs/{sample}.sort_bam.log"
  shell:
    "samtools sort {input} -o ${output} --threads 8"

rule index_bam:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.bam.bai"
  log: "logs/{sample}.index_bam.log"
  shell:
    "samtools index {input} {output}"

rule fix_mate_pairs: 
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.bam"    
  log: "logs/{sample}.fix_mate_pairs.log" 
  shell:  
    "java -jar $PICARD FixMateInformation\ I = {input} \ O = {output} \ SO = coordinate \ VALIDATION_STRINGENCY = LENINENT \ CREATE_INDEX = true"

rule filter_mapped_and_paired_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam"    
  log: "logs/{sample}.filter_mapped_and_paired_reads.log"
  shell:
    "bamtools filter \ -isMapped true \ -in {input} \ -out {output}"

rule remove_duplicate_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam"
  log: "logs/{sample}.remove_duplicate_reads.log"
  shell:
    "java -jar $PICARD MarkDuplicates \ INPUT = {input} \ OUTPUT = {output} \ VALIDATION_STRINGENCY=SILENT \ REMOVE_DUPLICATES = true \ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 4000"

rule add_read_groups:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.bam"
  log: "logs/{sample}.add_read_groups.log"
  shell:
    "java - jar $PICARD AddOrReplaceReadGroups \ INPUT = {input} \ OUTPUT = {output} \ RGLB = {wildcard.sample}.PE \ RGPL = Illumina \ RGPU = Group1 \ RGSM = {wildcard.sample}.PE" 

rule quality_filter_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.passed.bam"
  log: "logs/{sample}.quality_filter_reads.log"
  shell:
    "bamtools -filter - mapquality >=20 - in {input} - out {output}" 
