#!/usr/bin/env/ python
shell.prefix("set -o pipefail; ")
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
configfile: "config.yaml"

rule all:
    input:
        #trim reads
        expand("data/{sample}_{read}.fastq.gz", sample=config["SAMPLES"], read=["R1", "R2"]),
        #bwa_map_and_samtools_make_bam
        expand("data/trimmed_reads/{sample}_{read}_paired.fastq.gz", sample=config["SAMPLES"], read=["R1", "R2"]),
        #sort_bam
        expand("data/aligned_reads/{sample}.PE.bwa.bam", sample=config["SAMPLES"]),
        #index_bam
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.bam", sample=config["SAMPLES"]),
        #fix_mate_pairs
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.bam", sample=config["SAMPLES"]),
        #filter_mapped_and_paired_reads
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.bam", sample=config["SAMPLES"]),
        #remove_duplicate_reads
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam", sample=config["SAMPLES"]),
        #add_read_groups
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam", sample=config["SAMPLES"]),
        #quality_filter_reads
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.bam", sample=config["SAMPLES"]),
        #remove_louse_mitochondrial_reads
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.bam", sample=config["SAMPLES"]),
        #local_realignment_ID
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam", sample=config["SAMPLES"]),
        #local_realignment
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam", sample=config["SAMPLES"]),
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam.list", sample=config["SAMPLES"]),
        #call_variants
        expand("data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.bam", sample=config["SAMPLES"]),
        #combine_variant_files
        "data/variants/variants.list",
        #joint_variant_calling
        "data/variants/cohort.g.vcf.gz",
        #extract_snps
        "data/variants/combined.vcf.gz",
        #extract_indels
        "data/variants/combined.vcf.gz",
        #hard_filter_snps
        "data/variants/combined_snps.vcf",
        #hard_filter_indels
        "data/variants/combined_indels.vcf",
        #filter_maf
        "data/variants/combined_snps.vcf",
        #format_vcf_for_map_ped
        "data/variants/combined_snps.flt.maf.vcf",
        #make_map_ped
        "data/variants/combined_snps.flt.maf.no_chr.vcf",
        #filter_ld
        "data/variants/combined_snps.flt.maf.no_chr.plink",
        #convert_vcf_to_GESTE
        "data/variants/combined_snps.flt.maf.ld",
        "data/population_definitions.spid",
        #remove_temporary_files
        "data/variants/combined_snps.flt.maf.ld",
        "data/variants/combined_snps.flt.maf.no_chr.plink",
        #final
        "data/variants/combined_snps.flt.maf.ld.GESTE.txt"

rule trim_reads:
    input:
        expand("data/{{sample}}_{read}.fastq.gz", read=["R1", "R2"])
    output:
        "data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        "data/trimmed_reads/{sample}_R1_unpaired.fastq.gz",
        "data/trimmed_reads/{sample}_R2_paired.fastq.gz",
        "data/trimmed_reads/{sample}_R2_unpaired.fastq.gz"
    log: "logs/trim/{sample}.trim.log"
    shell:
        "module load trimmomatic;"
        "java -jar {config[TRIMMOMATIC]} PE -phred33 -trimlog {log} {input} {output} ILLUMINACLIP:{config[ADAPTERS]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    
rule bwa_map_and_samtools_make_bam:
    input:
        expand("data/trimmed_reads/{{sample}}_{read}_paired.fastq.gz", read=["R1", "R2"])
    output:
        "data/aligned_reads/{sample}.PE.bwa.bam"
    log: "logs/{sample}.map_and_bam.log"
    shell:
        "module load bwa samtools;"
        "bwa mem -t 8 {config[GENOME]} {input} | samtools view -Sb > {output}"

rule sort_bam:
    input:
        "data/aligned_reads/{sample}.PE.bwa.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.bam"
    log: "logs/{sample}.sort_bam.log"
    shell:
        "module load samtools;"
        "samtools sort {input} -o {output} --threads 8"

rule index_bam:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.bam.bai"
    log: "logs/{sample}.index_bam.log"
    shell:
        "module load samtools;"
        "samtools index {input} {output}"

rule fix_mate_pairs: 
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.bam"    
    log: "logs/{sample}.fix_mate_pairs.log" 
    shell: 
        "module load picard;"
        "java -jar {config[PICARD]} FixMateInformation INPUT={input} OUTPUT={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"

rule filter_mapped_and_paired_reads:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam"    
    log: "logs/{sample}.filter_mapped_and_paired_reads.log"
    shell:
        "module load bamtools;"
        "bamtools filter -isMapped true -in {input} -out {output}"

rule remove_duplicate_reads:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam"
    log: "logs/{sample}.remove_duplicate_reads.log"
    shell:
        "module load picard;"
        "java -jar {config[PICARD]} MarkDuplicates INPUT={input} OUTPUT={output} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE={log}"

rule add_read_groups:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.bam"
    log: "logs/{sample}.add_read_groups.log"
    shell:
        "module load picard;"
        "java -jar {config[PICARD]} AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={wildcards.sample}.PE RGPL=Illumina RGPU=Group1 RGSM={wildcards.sample}.PE" 

rule quality_filter_reads:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.bam"
    log: "logs/{sample}.quality_filter_reads.log"
    shell:
        "module load bamtools;"
        "bamtools filter -mapQuality '>=20' -in {input} -out {output}"

rule remove_louse_mitochondrial_reads:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam"
    log: "logs/{sample}.louse_mito_removal.log"
    shell:
        "module load samtools;"
        "samtools idxstats {input} | cut -f 1 | grep -v FJ* | xargs samtools view -b {input} > {output};"
        "samtools index -b {output} {output}.bai"

rule local_realignment_ID:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam.list"
    log: "logs/{sample}.local_realignment_ID.log"
    shell:
        "module load gatk;"
        "java -jar -Xmx4g {config[GATK]} -T RealignerTargetCreator --filter_mismatching_base_and_quals --num_threads 8 -R {config[GENOME]} -I {input} -o {output} --log_to_file {log}"

rule local_realignment:
    input:
        sample="data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam",
        list="data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.bam.list"
    output:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.bam"
    log: "logs/{sample}.local_realignment.log"
    shell:
        "module load gatk;"
        "java -jar -Xmx24g {config[GATK]} -I {input.sample} -R {config[GENOME]} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals {input.list} -o {output} --log_to_file {log}"

rule call_variants:
    input:
        "data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.bam"
    output:
        gvcf="data/variants/{sample}.g.vcf.gz",
        global_realign="data/aligned_reads/{sample}.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.global_realign.bam"
    log: "logs/{sample}.call_variants.log"
    shell:
        "module load gatk;"
        "java -jar -Xmx4g {config[GATK]} -T HaplotypeCaller -R {config[GENOME]} -I {input} -o {output.gvcf} -ERC GVCF -bamout {output.global_realign} --log_to_file {log}"

rule combine_variant_files:
    input:
        "data/variants/variants.list"
    output:
        "data/variants/cohort.g.vcf.gz"
    log: "logs/combine_variant_files.log"
    shell:
        "module load gatk;"
        "java -jar -Xmx4g {config[GATK]} -T CombineGVCFs -R {config[GENOME]} -V {input} -o {output} --log_to_file {log}"

rule joint_variant_calling:
    input:
        "data/variants/cohort.g.vcf.gz"
    output:
        "data/variants/combined.vcf.gz"
    log: "logs/joint_variant_calling.log"
    shell:
        "module load gatk;"
        "java -jar -Xmx4g {config[GATK]} -T GenotypeGVCFs -R {config[GENOME]} -V {input} -o {output} --log_to_file {log}"

rule extract_snps:
    input:
        "data/variants/combined.vcf.gz"
    output:
        "data/variants/combined_snps.vcf"
    log: "logs/extract_snps.log"
    shell:
        "module load gatk;"
        "java -jar {config[GATK]} -T SelectVariants -R {config[GENOME]} -V {input} -o {output} -selectType SNP"

rule extract_indels:
    input:
        "data/variants/combined.vcf.gz"
    output:
        "data/variants/combined_indels.vcf"
    log: "logs/extract_indels.log"
    shell:
        "module load gatk;"
        "java -jar {config[GATK]} -T SelectVariants -R {config[GENOME]} -V {input} -o {output} -selectType INDEL"

rule hard_filter_snps:
    input: 
        "data/variants/combined_snps.vcf" 
    output:
        "data/variants/combined_snps.flt.vcf"
    log: "logs/hard_filter_snps.log"
    shell:
        "module load gatk;"
        "java -jar {config[GATK]} -T VariantFiltration -R {config[GENOME]} -V {input} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || HaplotypeScore > 13.0' --filterName 'NGS_snp_filter' -o {output} --missingValuesInExpressionsShouldEvaluateAsFailing"

rule hard_filter_indels:
    input:
        "data/variants/combined_indels.vcf"
    output:
        "data/variants/combined_indels.vcf"
    log: "logs/hard_filter_indels.log"
    shell:
        "module load gatk;"
        "java -jar {config[GATK]} -T VariantFiltration -R {config[GENOME]} -V {input} --filterExpression 'QD < 2.0 || FS > 200.0 || MQ < 40.0 || ReadPosRankSum < -20.0' --filterName 'NGS_indel_filter' -o {output} --missingValuesInExpressionsShouldEvaluateAsFailing"

rule filter_maf:
    input:
        "data/variants/combined_snps.flt.vcf"
    output:
        "data/variants/combined_snps.flt.maf.vcf"
    log: "logs/filter_maf.log"
    shell:
        "module load gatk;"
        "java -jar {config[GATK]} -T SelectVariants -R {config[GENOME]} -V {input} -o {output} -select 'AF > 0.05' --restrictAllelesTo BIALLELIC"

rule format_vcf_for_map_ped:
    input:
        "data/variants/combined_snps.flt.maf.vcf"
    output:
        "data/variants/combined_snps.flt.maf.no_chr.vcf"
    log: "logs/format_vcf_for_map_ped.log"
    shell:
        "sed 's/^DS//' {input} | sed '/^FJ/d' | sed 's/ID=DS/ID=/' | sed '/ID=FJ/d'  > {output}"

rule make_map_ped:
    input:
        "data/variants/combined_snps.flt.maf.no_chr.vcf"
    output:
        "data/variants/combined_snps.flt.maf.no_chr.plink"
    output:
    log: "logs/make_bed_map.log"
    shell:
        "module load vcftools;"
        "{config[VCFTOOLS]} --vcf {input} --plink --out {output};"
        " > {output}"

rule filter_ld:
    input:
        "data/variants/combined_snps.flt.maf.no_chr.plink"
    output:
        "data/variants/combined_snps.flt.maf.ld"
    log: "logs/filter_ld.log"
    shell:
        "module load plink;"
        "plink --file {input} --allow-extra-chr --indep 50 5 1.25;"
        "plink --file {input} --allow-extra-chr --extract plink.prune.in --make-bed;"
        "plink --bfile plink --allow-extra-chr --recode vcf --out {output};"
        " > {output};"
        "rm plink*"

rule convert_vcf_to_GESTE:
    input:
        vcf="data/variants/combined_snps.flt.maf.ld",
        spid="data/population_definitions.spid"
    output:
        "data/variants/combined_snps.flt.maf.ld.GESTE.txt"
    log: "logs/convert_vcf_to_GESTE.log"
    shell:
        "module load pdgspider;"
        "java -Xmx1024m -Xms512M -jar {config[PDGSPIDER]} -inputfile {input.vcf}.vcf -inputformat VCF -outputfile {output} -outputformat GESTE_BAYE_SCAN -spid {input.spid}"

rule remove_temporary_files:
    input:
        "data/variants/combined_snps.flt.maf.ld",
        "data/variants/combined_snps.flt.maf.no_chr.plink"
    log: "logs/remove_temporary_files.log"
    shell:
        "rm {input}"