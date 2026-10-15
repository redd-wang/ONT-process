#!/bin/bash

#SBATCH --job-name=ONT-methylation-call
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=methylation_call_5mC.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

##this script starts with uBAM directly from Dorado, unaligned to ref genome and EM-seq cov files##

echo "loading modkit...calling 5mC"
~/peteslab/Redd/bin/modkit/modkit pileup align/FGC.bam \
FGC_nanopore_5mC.bed \
--modified-bases 5mC \
--ref /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
--log 5mC_modkit_log.txt \
--cpg --threads 16

echo "loading modkit...calling 5hmC"
~/peteslab/Redd/bin/modkit/modkit pileup align/FGC.bam \
FGC_nanopore_5hmC.bed \
--modified-bases 5hmC \
--ref /shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
--log 5hmC_modkit_log.txt \
--cpg --threads 16