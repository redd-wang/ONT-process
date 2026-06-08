#!/bin/bash

#SBATCH --job-name=assembly-35CR
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --output=process.out
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

##This script is for ONT cleanup and assembly using fastq which is converted from uBAM, output from Dorato##

echo "de novo assembly by hifiasm, only ONT"

~/peteslab/Redd/bin/hifiasm-0.25.0/hifiasm -o 35CR_run1.asm --ont -t32 35CR_run1.fastq
