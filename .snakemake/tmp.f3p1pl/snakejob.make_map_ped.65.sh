#!/bin/sh
# properties = {"jobid": 65, "rule": "make_map_ped", "output": ["data/variants/combined_snps.maf.no_chr.plink"], "threads": 1, "input": ["data/variants/combined_snps.maf.no_chr.vcf"], "resources": {}, "local": false, "cluster": {}, "params": {}}
cd /gpfs/scratch/cah422/NGS_map-snakemake && \
/opt/rh/python33/root/usr/bin/python3 -m snakemake data/variants/combined_snps.maf.no_chr.plink --snakefile /gpfs/scratch/cah422/NGS_map-snakemake/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files data/variants/combined_snps.maf.no_chr.vcf /gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.f3p1pl --latency-wait 5 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules make_map_ped  && touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.f3p1pl/65.jobfinished" || (touch "/gpfs/scratch/cah422/NGS_map-snakemake/.snakemake/tmp.f3p1pl/65.jobfailed"; exit 1)

