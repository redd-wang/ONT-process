#!/bin/bash

#SBATCH --job-name=FGC-fastq-split
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=FGC-T2T-alignment.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

#separate one fastq to multiple files for parallel processing
module load SeqKit
echo "partition fastq to 10 segments"
mkdir -p split
seqkit split2 -p 10 \
-j 16 \
-O split_fastq \
-e .gz \
fastq/FGC.fastq