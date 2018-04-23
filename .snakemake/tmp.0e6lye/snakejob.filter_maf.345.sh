#!/bin/sh
# properties = {"params": {}, "resources": {}, "local": false, "output": ["data/variants/combined_snps.maf.vcf"], "input": ["data/variants/combined_snps.vcf"], "threads": 1, "jobid": 345, "cluster": {}, "rule": "filter_maf"}
cd /gpfs/scratch/cah422/NGS_map-snakemake && \
/opt/rh/python33/root/usr/bin/python3 -m snakemake data/variants/combined_snps.maf.vcf --snakefile /gpfs/scratch/cah422/NGS_map-snakemake/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files data/variants/combined_snps.vcf /gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.0e6lye --latency-wait 5 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules filter_maf  && touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.0e6lye/345.jobfinished" || (touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.0e6lye/345.jobfailed"; exit 1)

