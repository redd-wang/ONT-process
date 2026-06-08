#!/bin/bash

#SBATCH --job-name=CNV-calling-35CR
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=CNV-calling.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

echo "call CNV in the scaffold-aligned bam files 35CR with correct chromosome number"
module load CNVkit

cnvkit.py batch \
    align/Lu35CR_aligned.bam \
    --seq-method wgs \
    --drop-low-coverage \
    --reference /home/rwang4/CoRAL/reference/hg38full_ref_5k.cnn \
    --scatter --diagram \
    -d cnv_output/ \
    -p 8