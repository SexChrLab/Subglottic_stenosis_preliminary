#!/bin/bash
#SBATCH --job-name=rna  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH -q tempboost

snakemake --snakefile process_rna.snakefile -j 19 --cluster "sbatch --mem=48000 -q tempboost -c 4 -t 40:00:00" --rerun-incomplete
