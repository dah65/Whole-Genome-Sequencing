#!/bin/sh
# properties = {"output": ["data/variants/combined_snps.maf.ld.vcf"], "cluster": {}, "jobid": 103, "resources": {}, "input": ["data/variants/combined_snps.maf.vcf"], "local": false, "threads": 1, "params": {}, "rule": "filter_ld"}
cd /gpfs/scratch/cah422/NGS_map-snakemake && \
/opt/rh/python33/root/usr/bin/python3 -m snakemake data/variants/combined_snps.maf.ld.vcf --snakefile /gpfs/scratch/cah422/NGS_map-snakemake/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files data/variants/combined_snps.maf.vcf /gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.fod0yl --latency-wait 5 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules filter_ld  && touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.fod0yl/103.jobfinished" || (touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.fod0yl/103.jobfailed"; exit 1)

