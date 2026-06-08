#!/bin/bash
#SBATCH --job-name=samtools_darkregion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=darkregiongenes_call.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END
#SBATCH --output=logs/05_depth.log

set -euo pipefail

module load SAMtools

mkdir -p depth logs

echo "[$(date)] Running samtools depth for Lu35CR_T2T"

# ── Per-region depth over dark gene BED ─────────────────────
#samtools depth \
 #       -b /home/rwang4/peteslab/Redd/reference/chm13/dark_region_genes.bed \
 #       -q 1 \
 #       --threads 8 \
 #       alignT2T/35CR_T2Taligned.bam > depth/Lu35CR_T2Taligned_depth.tsv

# Genome-wide depth for normalisation
samtools coverage \
    --min-MQ 1 \
    alignT2T/35CR_T2Taligned.bam \
    > depth/Lu35CR_genome.cov

module load fhPython
python3 05_coverage_samtools_depth.py \
    --depth  depth/Lu35CR_T2Taligned_depth.tsv \
    --bed    /home/rwang4/peteslab/Redd/reference/chm13/dark_region_genes.bed \
    --out    depth/Lu35CR_dark.regions.bed \
    --genome_cov depth/Lu35CR_genome.cov