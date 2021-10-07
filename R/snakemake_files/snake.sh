#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1440
#SBATCH --mem=1000

module load python
snakemake -j 100 --latency-wait 30 --cluster "sbatch --mem=10000 -N 1 -n 4"
