#!/bin/bash

#SBATCH --job-name=ONT-SV-calling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --output=SV_call.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

module load fhPython
module load Sniffles/2.6.2-gfbf-2024a

##calling somatic SVs
#echo "calling somatic SVs"
#sniffles --input ../35CR_run1_aligned.bam \
#--vcf 35CR_R1_somaticSV.vcf \
#--threads 16 \
#--mosaic

##calling germline SVs
echo "calling germline SVs"
sniffles --input alignT2T/35CR_T2Taligned.bam \
--vcf SV/35CR_T2T_germlineSV.vcf \
--reference /home/rwang4/peteslab/Redd/reference/chm13/chm13v2.0.fa \
--threads 16 

