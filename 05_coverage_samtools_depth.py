"""
coverage_samtools_depth.py
--------------------------
Summarizes samtools depth per-base output into per-gene mean coverage,
replacing mosdepth regions output.

Produces the same format as mosdepth .regions.bed.gz so the downstream
06_intersect_annotate.py script works unchanged.

Usage:
    # Step 1: run samtools depth on cluster
    samtools depth \
        -b dark_region_genes.bed \
        -q 1 \
        --threads 8 \
        your.bam > sample_depth.tsv

    # Step 2: summarize to per-gene
    python coverage_samtools_depth.py \
        --depth  sample_depth.tsv \
        --bed    dark_region_genes.bed \
        --out    sample_dark.regions.bed

    # Run for both normal and PDX:
    python coverage_samtools_depth.py --depth normal_depth.tsv --bed dark_region_genes.bed --out mosdepth_out/normal_dark.regions.bed
    python coverage_samtools_depth.py --depth PDX_depth.tsv    --bed dark_region_genes.bed --out mosdepth_out/PDX_dark.regions.bed

Dependencies:
    pip install --user pandas numpy
"""

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict


def load_bed(bed_path):
    """Load BED file into DataFrame."""
    genes = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue
            genes.append({
                "chrom": cols[0],
                "start": int(cols[1]),
                "end":   int(cols[2]),
                "name":  cols[3],
            })
    return pd.DataFrame(genes)


def summarize_depth(depth_tsv, bed_df, out_path):
    """
    Read samtools depth output (chrom, pos, depth per line)
    and compute mean depth per gene region.

    samtools depth output format:
        chr1    1000001    25
        chr1    1000002    27
        ...
    Positions with 0 depth are NOT reported by default —
    we account for this by using gene length as denominator.
    """

    print("Loading BED regions...")
    # Build interval lookup: chrom -> sorted list of (start, end, name)
    intervals = defaultdict(list)
    for _, row in bed_df.iterrows():
        intervals[row["chrom"]].append((row["start"], row["end"], row["name"]))

    # Accumulators: gene_name -> (sum_depth, gene_length)
    gene_sum    = defaultdict(int)
    gene_length = {row["name"]: (row["end"] - row["start"])
                   for _, row in bed_df.iterrows()}
    gene_coords = {row["name"]: (row["chrom"], row["start"], row["end"])
                   for _, row in bed_df.iterrows()}

    print("Reading samtools depth output...")
    chunk_size = 500_000
    chunks_read = 0
    current_chrom = None
    chrom_intervals = []

    with open(depth_tsv) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            pos   = int(cols[1])   # 1-based in samtools depth
            depth = int(cols[2])

            # Update interval list when chromosome changes
            if chrom != current_chrom:
                current_chrom  = chrom
                chrom_intervals = sorted(intervals.get(chrom, []))

            # Find which gene(s) this position falls in
            # Simple linear scan — efficient enough for ONT depth files
            # (genes are non-overlapping in BED)
            for (gstart, gend, gname) in chrom_intervals:
                if pos - 1 >= gstart and pos - 1 < gend:   # convert to 0-based
                    gene_sum[gname] += depth
                elif gstart > pos:
                    break   # sorted, no need to continue

            chunks_read += 1
            if chunks_read % 5_000_000 == 0:
                print(f"  Processed {chunks_read:,} positions...")

    print(f"  Total positions processed: {chunks_read:,}")

    # Compute mean depth per gene
    # Positions with 0 depth are absent from samtools depth output
    # → use gene length as denominator (counts zero positions correctly)
    print("Computing per-gene mean depth...")
    results = []
    for _, row in bed_df.iterrows():
        gname  = row["name"]
        glen   = gene_length[gname]
        gsum   = gene_sum.get(gname, 0)
        mean_d = gsum / glen if glen > 0 else 0.0
        results.append({
            "chrom": row["chrom"],
            "start": row["start"],
            "end":   row["end"],
            "name":  gname,
            "mean_depth": round(mean_d, 3),
        })

    out_df = pd.DataFrame(results)

    # Save in mosdepth-compatible format: chrom start end name depth
    out_df.to_csv(out_path, sep="\t", header=False, index=False)
    print(f"Saved: {out_path}  ({len(out_df)} genes)")

    # Quick summary
    print(f"\nCoverage summary:")
    print(f"  Median depth:         {out_df['mean_depth'].median():.1f}x")
    print(f"  Mean depth:           {out_df['mean_depth'].mean():.1f}x")
    print(f"  Genes with 0x depth:  {(out_df['mean_depth'] == 0).sum()}")
    print(f"  Genes with <5x depth: {(out_df['mean_depth'] < 5).sum()}")
    print(f"  Genes with >30x:      {(out_df['mean_depth'] > 30).sum()}")

    return out_df


def compute_genome_median(bam_path, sample_name, out_path, threads=8):
    """
    Print the samtools command to get genome-wide median depth.
    We use samtools coverage (fast, gives mean per chromosome)
    then take the weighted median.
    """
    print(f"\nTo get genome-wide median depth for normalisation, run:")
    print(f"  samtools coverage --min-MQ 1 {bam_path} > {out_path}")
    print(f"  Then pass --genome_cov {out_path} to this script")


def load_genome_coverage(genome_cov_path):
    """
    Parse samtools coverage output to get genome-wide median depth.
    samtools coverage columns:
    #rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    """
    rows = []
    with open(genome_cov_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 7:
                continue
            chrom    = cols[0]
            length   = int(cols[2]) - int(cols[1]) if cols[2].isdigit() else 0
            meandepth = float(cols[6])
            rows.append({"chrom": chrom, "length": length, "meandepth": meandepth})
    df = pd.DataFrame(rows)
    # Weighted median across chromosomes
    if len(df) == 0:
        return 30.0   # fallback
    total_len = df["length"].sum()
    df = df.sort_values("meandepth")
    df["cumlen"] = df["length"].cumsum()
    median_row = df[df["cumlen"] >= total_len / 2].iloc[0]
    return float(median_row["meandepth"])


def main():
    parser = argparse.ArgumentParser(
        description="Summarize samtools depth to per-gene coverage (mosdepth replacement)")
    parser.add_argument("--depth",      required=True,
                        help="samtools depth output TSV")
    parser.add_argument("--bed",        required=True,
                        help="Dark region gene BED file")
    parser.add_argument("--out",        required=True,
                        help="Output per-gene depth BED")
    parser.add_argument("--genome_cov", default=None,
                        help="samtools coverage output for genome median (optional)")
    args = parser.parse_args()

    bed_df = load_bed(args.bed)
    print(f"BED regions loaded: {len(bed_df)} genes")

    out_df = summarize_depth(args.depth, bed_df, args.out)

    if args.genome_cov:
        median = load_genome_coverage(args.genome_cov)
        print(f"\nGenome-wide median depth: {median:.1f}x")
        print(f"Ratio (gene median / genome): "
              f"{out_df['mean_depth'].median() / median:.2f}x")


if __name__ == "__main__":
    main()
