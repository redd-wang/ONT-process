#!/bin/bash

#SBATCH --job-name=dark-region-annotation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=darkregionannotation.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

module load BEDTools
module load fhPython

python3 06_intersect_annotate2.py \
        --dark_bed      novel_t2t_genes.bed \
        --pdx_vcf       SV/35CR_T2T_germlineSV.vcf \
        --pdx_cov       depth/Lu35CR_dark.regions.bed \
        --pdx_genome    depth/Lu35CR_genome.cov \
        --out_report    dark_genes_report.tsv \
        --out_hits      dark_genes_hits.tsv


python3 07_plot_dark_region_results2.py \
    --report dark_genes_report.tsv \
    --out    dark_region_novel_gene_results.png

