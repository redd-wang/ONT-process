#!/bin/bash

#SBATCH --job-name=ONT-align-T2T
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=T2T-CHM13-align.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

#align uBAM to ref genome
module load minimap2
module load SAMtools

mkdir -p alignT2T
echo "aligning fastq to the T2T-CHM13 genome with methylation tags"
minimap2 -ax map-ont -y -t 16 \
/home/rwang4/peteslab/Redd/reference/chm13v2.0.fa \
fastq/35CR_merge.fastq | \
samtools sort -@ 16 -o alignT2T/35CR_T2Taligned.bam
samtools index alignT2T/35CR_T2Taligned.bam

