#!/bin/bash

#SBATCH --job-name=ONT-FGC-hg38align
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=FGC_hg38_alignment.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

##this script starts with uBAM directly from Dorado##

#align uBAM to ref genome
module load minimap2
module load SAMtools

echo "merge all sub run files to single bam"
samtools cat -o unalignedBAM/FGC.bam unalignedBAM/*.bam
samtools view -H unalignedBAM/FGC.bam | grep '@RG'
samtools flagstat unalignedBAM/FGC.bam

echo "converting uBAM file from Dorado to fastq with MM/ML tags"
samtools fastq -T MM,ML \
unalignedBAM/FGC.bam \
> fastq/FGC.fastq

echo "aligning fastq to the hg38 genome with methylation tags"
minimap2 -ax map-ont -y -t 16 \
/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
fastq/FGC.fastq | \
samtools sort -@ 16 -o align/FGC.bam
samtools index align/FGC.bam


