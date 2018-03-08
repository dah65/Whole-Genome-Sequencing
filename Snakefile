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
    expand("genomes/{genome}", genome=config["GENOME"]),
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
    """
    java -jar $PICARD FixMateInformation  
      I = {input} 
      O = {output}  
      SO = coordinate  
      VALIDATION_STRINGENCY = LENINENT  
      CREATE_INDEX = true"
    """

rule filter_mapped_and_paired_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam"    
  log: "logs/{sample}.filter_mapped_and_paired_reads.log"
  shell:
    """
    bamtools filter  
      -isMapped true  
      -in {input} 
      -out {output}
    """

rule remove_duplicate_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam"
  log: "logs/{sample}.remove_duplicate_reads.log"
  shell:
    """
    java -jar $PICARD MarkDuplicates  
      INPUT = {input} 
      OUTPUT = {output}  
      VALIDATION_STRINGENCY=SILENT  
      REMOVE_DUPLICATES = true  
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 4000"
    """

rule add_read_groups:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.bam"
  log: "logs/{sample}.add_read_groups.log"
  shell:
    """
    java - jar $PICARD AddOrReplaceReadGroups  
      INPUT = {input} 
      OUTPUT = {output}  
      RGLB = {wildcard.sample}.PE  
      RGPL = Illumina  
      RGPU = Group1 
      RGSM = {wildcard.sample}.PE 
    """

rule quality_filter_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.passed.bam"
  log: "logs/{sample}.quality_filter_reads.log"
  shell:
    "bamtools -filter - mapquality >=20 - in {input} - out {output}"

rule remove_louse_mitochondrial_reads:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.bam"
  log: "logs/{sample}.louse_mito_removal.log"
  shell:
    "samtools idxstats {input} | cut -f 1 | grep -v FJ* | xargs samtools view -b {input} > {output}"

rule local_realignment_ID:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.bam", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.bam.list"
  log: "logs/{sample}.local_realignment_ID.log"
  shell:
    """
    java - jar - Xmx4g -jar $GATK  
      -T RealignerTargetCreator 
      --filter_mismatching_base_and_quals  
      --num_threads 8  
      -I {input}  
      -o {output}
    """

rule local_realignment:
  input:
    sample=expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.bam", sample=config["SAMPLES"]),
    list=expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.bam.list", sample=config["SAMPLES"])
  output:
    "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.realign.bam"
  log: "logs/{sample}.local_realignment.log"
  shell:
    """
    java -jar -Xmx24g $GATK  
      -I {input.sample}  
      -R {config[GENOME]}  
      --filter_mismatching_base_and_quals  
      -T IndelRealigner  
      -targetIntervals {input.list}  
      -o {output}
    """

rule call_snps:
  input:
    expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.mito_removed.realign.bam", sample=config["SAMPLES"])
  output:
    "data/variants/{sample}.bcf"
  log: "logs/{sample}.call_snps.log"
  shell:
    """
    samtools mpilup 
      -C50 
      -ugf {config[GENOME]}  
      {input} |  
      bcftools call -vmO b > {output}
    """

rule filter_snps:
  input:
    expand("data/variants/{sample}.bcf", sample=config["SAMPLES"])
  output:
    "data/variants/{sample}.flt.vcf"
  log: "logs/{sample}.filter_snps.log"
  shell:
    "bcftools veiw {input} | vcfutils.pl varFilter -d 5 -D 100 > {output}"

rule index_and_compress_vcf:
  input:
    expand("data/variants/{sample}.flt.vcf", sample=config["SAMPLES"])
  output:
    "data/variants/{sample}.flt.vcf.gz"
  log: "logs/{sample}.index_and_compress_vcf.log"
  run:
    "bgzip -c {input} | tabix -p vcf {output}"

rule get_snp_stats:
  input:
    expand("data/variants/{sample}.flt.vcf.gz", sample=config["SAMPLES"])
  output:
    "data/variants/{sample}.flt.vcf.stats.txt"
  log: "logs/{sample}.get_snp_stats.log"
  shell:
    "vcf-stats {input} > {output}"

rule call_consensus:
  input:
    expand("data/variants/{sample}.flt.vcf", sample=config["SAMPLES"])
  output:
    "data/consensus_sequences/{sample}.bwa.consensus.fq.gz"
  log: "logs/{sample}.call_consensus.log"
  shell:
    "vcfutils.pl vcf2fq -d 5 -D 100 {input} | gzip > {output}"
