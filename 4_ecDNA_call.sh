#!/bin/bash

#SBATCH --job-name=ONT-ecDNA-calling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=ecDNA_call.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

##load conda base environment##
source ~/.bashrc
##load conda coral environment##
conda activate coral_env

##set paths
BAM="/home/rwang4/peteslab/Redd/ONT/Lu35CR/align/Lu35CR_aligned.bam"
CNS="/home/rwang4/peteslab/Redd/ONT/Lu35CR/cnv_output/Lu35CR_aligned.cns"
CENTRO="/home/rwang4/CoRAL/coral/supplemental_data/GRCh38_centromere.bed"

echo "Seed step starts"
coral seed \
--cn-seg $CNS \
--output-prefix ecDNA/Lu35CR \
--lr-bam $BAM \
--centromere-file $CENTRO \

echo "restruction step starts"
coral reconstruct \
--lr-bam $BAM \
--cnv-seed ecDNA/Lu35CR_CNV_SEEDS.bed \
--cn-seg $CNS \
--solver scip \
--output-prefix ecDNA/Lu35CR \
--solver-threads 16

##deactivate coral
conda deactivate
##deactivate base
conda deactivate
module purge


