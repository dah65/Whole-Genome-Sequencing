#!/bin/sh
# properties = {"input": ["data/aligned_reads/Col_Body_1.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.bam"], "output": ["data/aligned_reads/Col_Body_1.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.global_realign.bam", "data/variants/Col_Body_1.g.vcf.gz"], "local": false, "jobid": 317, "rule": "call_variants", "cluster": {}, "resources": {}, "params": {}, "threads": 1}
cd /gpfs/scratch/cah422/NGS_map-snakemake && \
/opt/rh/python33/root/usr/bin/python3 -m snakemake data/aligned_reads/Col_Body_1.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.global_realign.bam data/variants/Col_Body_1.g.vcf.gz --snakefile /gpfs/scratch/cah422/NGS_map-snakemake/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files data/aligned_reads/Col_Body_1.PE.bwa.sorted.fixed.filtered.postdup.RG.passed.mito_removed.local_realign.bam /gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.72gsrd --latency-wait 5 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules call_variants  && touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.72gsrd/317.jobfinished" || (touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.72gsrd/317.jobfailed"; exit 1)

