#!/bin/bash

#SBATCH --job-name=ONT-methylation-cleanup
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=methylation_call_5hmC.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

##this script starts with uBAM directly from Dorado, unaligned to ref genome and EM-seq cov files##

#align uBAM to ref genome
module load minimap2
module load SAMtools
echo "converting uBAM files from Dorado to fastq with MM/ML tags"
samtools fastq -T MM,ML \
../PBI57575_pass_3e87bdd2_8bdefc50_0.bam \
> 35CR_run1.fastq
echo "aligning fastq to the hg38 genome with methylation tags"
minimap2 -ax map-ont -y -t 16 \
/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
35CR_run1.fastq | \
samtools sort -@ 16 -o 35CR_run1_aligned.bam
samtools index 35CR_run1_aligned.bam

#methylation output to bed files
echo "loading modkit...calling 5mC"
~/peteslab/Redd/bin/modkit/modkit pileup 35CR_run1_aligned.bam \
35CR_run1_nanopore_5mC.bed \
--modified-bases 5mC \
--ref /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
--log 5mC_modkit_log.txt \
--cpg --threads 16

echo "loading modkit...calling 5hmC"
~/peteslab/Redd/bin/modkit/modkit pileup 35CR_run1_aligned.bam \
35CR_run1_nanopore_5hmC.bed \
--modified-bases 5hmC \
--ref /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
--log 5hmC_modkit_log.txt \
--cpg --threads 16