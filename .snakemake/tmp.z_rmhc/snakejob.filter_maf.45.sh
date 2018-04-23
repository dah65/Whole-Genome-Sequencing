#!/bin/sh
# properties = {"output": ["data/variants/combined_snps.maf.vcf"], "rule": "filter_maf", "cluster": {}, "jobid": 45, "local": false, "input": ["data/variants/combined_snps.vcf"], "resources": {}, "params": {}, "threads": 1}
cd /gpfs/scratch/cah422/NGS_map-snakemake && \
/opt/rh/python33/root/usr/bin/python3 -m snakemake data/variants/combined_snps.maf.vcf --snakefile /gpfs/scratch/cah422/NGS_map-snakemake/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files data/variants/combined_snps.vcf /gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.z_rmhc --latency-wait 5 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules filter_maf  && touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.z_rmhc/45.jobfinished" || (touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.z_rmhc/45.jobfailed"; exit 1)

