#!/bin/bash

#SBATCH --job-name=ONT-EMseq_correlation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --output=methylation_correlation.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

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

echo "calculate spearman R between ONT 5hmC calling with EM-seq of cell line"

module load fhPython

python3 <<EOF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

def run_validation():
    print("Loading data for plotting...")
    
    # 1. Load Nanopore
    # Column 9 = Coverage, Column 10 = Percent Methylated in most modkit versions
    nano = pd.read_csv("35CR_run1_nanopore_5hmC.bed", sep="\t", header=None, 
                       usecols=[0, 1, 9, 10], names=["chrom", "pos", "cov_nano", "freq_nano"])
    
    # Standardize chromosome names (strips 'chr' if present so both match)
    nano['chrom'] = nano['chrom'].astype(str).str.replace('chr', '')
    nano['freq_nano'] = nano['freq_nano'] / 100.0

    # 2. Load EM-seq
    em_file = "nsorted.human.35CR.HMC_AZE036_CL_EMSeq.bismark.cov.gz"
    emseq = pd.read_csv(em_file, sep="\t", header=None, compression='gzip',
                        usecols=[0, 1, 3], names=["chrom", "pos", "freq_em"])
    
    emseq['chrom'] = emseq['chrom'].astype(str).str.replace('chr', '')
    emseq['pos'] = emseq['pos'] - 1 # Bismark is 1-based, Modkit is 0-based
    emseq['freq_em'] = emseq['freq_em'] / 100.0

    # 3. Merge and filter
    combined = pd.merge(nano, emseq, on=["chrom", "pos"])
    df = combined[combined['cov_nano'] >= 10].dropna()

    if df.empty:

        print("Error: No overlapping sites found. Check if your reference names match!")
        return

    r_val, _ = spearmanr(df['freq_nano'], df['freq_em'])
    
    # --- PLOTTING ---
    plt.figure(figsize=(8, 8))
    hb = plt.hexbin(df['freq_em'], df['freq_nano'], gridsize=50, cmap='YlGnBu', bins='log')
    plt.colorbar(hb, label='log10(count of CpGs)')
    plt.plot([0, 1], [0, 1], color='red', linestyle='--', alpha=0.5, label='Identity')
    
    plt.title(f'Methylation Correlation: Nanopore vs EM-seq\n(Spearman R = {r_val:.3f})')
    plt.xlabel('EM-seq_CL (Bismark)')
    plt.ylabel('Nanopore 5hmC (Modkit)')
    plt.legend()
    
    plt.savefig("35CR_run1_spearman_plot_5hmC_CL_cutoff10.pdf")
    print("Plot saved to 35CR_run1_spearman_plot_5hmC_CL_cutoff10.pdf")

run_validation()

EOF