#!/bin/bash

#SBATCH --job-name=ONT-EMseq_correlation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --output=methylation_correlation.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

set -euo pipefail

echo "calculate spearman R between ONT 5mC calling with EM-seq of cell line"

module load fhPython

python3 <<EOF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

def run_validation():
    print("Loading data for plotting...")
    
    # 1. Load Nanopore (removed $ from filenames)
    # Column 9 = Coverage, Column 10 = Percent Methylated in most modkit versions
    nano = pd.read_csv("35CR_run1_nanopore.bed", sep="\t", header=None, 
                       usecols=[0, 1, 9, 10], names=["chrom", "pos", "cov_nano", "freq_nano"])
    
    # Standardize chromosome names (strips 'chr' if present so both match)
    nano['chrom'] = nano['chrom'].astype(str).str.replace('chr', '')
    nano['freq_nano'] = nano['freq_nano'] / 100.0

    # 2. Load EM-seq
    em_file = "nsorted.human.35CR.EM_AZE036_CL_EMSeq.bismark.cov.gz"
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
    plt.ylabel('Nanopore (Modkit)')
    plt.legend()
    
    plt.savefig("35CR_run1_spearman_plot_10cutoff.pdf")
    print("Plot saved to 35CR_run1_correlation_plot_10cutoff.pdf")

run_validation()

EOF