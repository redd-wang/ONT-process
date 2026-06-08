"""
07_plot_dark_region_results.py
------------------------------
Visualizes dark region gene alteration findings:
  Panel 1 — Genome-wide coverage ratio (PDX vs normal) highlighting dark genes
  Panel 2 — Altered dark gene heatmap (gene × evidence type)
  Panel 3 — Coverage ratio scatter: normal vs PDX per dark gene
  Panel 4 — SV type breakdown for SVs overlapping dark genes

Usage:
    python 07_plot_dark_region_results.py \
        --report dark_genes_report.tsv \
        --out    dark_region_results.png

Dependencies:
    pip install --user pandas matplotlib seaborn numpy
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

SV_COLORS = {
    "DEL": "#E63946", "INS": "#457B9D",
    "BND": "#9B2226", "DUP": "#2A9D8F", "INV": "#E9C46A",
}


def load(report_path):
    df = pd.read_csv(report_path, sep="\t")
    for col in ["normal_cov_ratio","pdx_cov_ratio","normal_depth","pdx_depth"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def plot_results(df, out):
    fig = plt.figure(figsize=(18, 16))
    fig.suptitle("Dark Region Gene Alterations — PDX vs Normal\n(T2T-CHM13 mapped)",
                 fontsize=14, fontweight="bold", y=0.99)

    gs = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, :])    # genome-wide coverage — full width
    ax2 = fig.add_subplot(gs[1, 0])    # heatmap
    ax3 = fig.add_subplot(gs[1, 1])    # scatter

    # ── Panel 1: Per-chromosome coverage ratio ─────────────
    ax = ax1
    ax.set_facecolor("#f8f9fa")
    ax.spines[["top","right"]].set_visible(False)

    chrom_order_present = [c for c in CHROM_ORDER if c in df["chrom"].values]
    x_pos = []
    x_labels = []
    cursor = 0
    chrom_gaps = 0.5

    for chrom in chrom_order_present:
        sub = df[df["chrom"] == chrom].copy()
        if len(sub) == 0:
            cursor += chrom_gaps
            continue

        sub["x"] = cursor + np.linspace(0, 0.8, len(sub))
        x_pos.append(cursor + 0.4)
        x_labels.append(chrom.replace("chr",""))

        for _, row in sub.iterrows():
            ratio = row["pdx_cov_ratio"]
            if np.isnan(ratio):
                continue
            color = "#E63946" if ratio < 0.5 else \
                    "#2A9D8F" if ratio > 2.5 else \
                    "#aaaaaa"
            size  = 40 if row["cancer_relevant"] else 15
            alpha = 0.9 if row["cancer_relevant"] else 0.5
            marker = "*" if row["cancer_relevant"] else "o"
            ax.scatter(row["x"], ratio, c=color, s=size,
                       alpha=alpha, marker=marker,
                       zorder=3 if row["cancer_relevant"] else 2)

            if row["cancer_relevant"] and (ratio < 0.5 or ratio > 2.5):
                ax.annotate(row["gene"],
                            xy=(row["x"], ratio),
                            xytext=(row["x"], ratio + 0.15),
                            fontsize=6.5, ha="center", color=color,
                            fontweight="bold",
                            arrowprops=dict(arrowstyle="-", lw=0.4, color=color))

        cursor += 1 + chrom_gaps

    ax.axhline(1.0, color="#333", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.axhline(0.5, color="#E63946", linewidth=0.6, linestyle=":", alpha=0.7)
    ax.axhline(2.5, color="#2A9D8F", linewidth=0.6, linestyle=":", alpha=0.7)
    ax.text(cursor, 0.52, "del threshold", fontsize=7, color="#E63946")
    ax.text(cursor, 2.52, "dup threshold", fontsize=7, color="#2A9D8F")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.set_ylabel("PDX Coverage Ratio\n(vs genome median)", fontsize=10)
    ax.set_title("Dark Region Gene Coverage — PDX/Normal Depth Ratio  (★ = cancer-relevant gene)",
                 fontsize=10, loc="left")
    ax.set_ylim(-0.1, max(4.0, df["pdx_cov_ratio"].dropna().quantile(0.99) + 0.5))

    red_patch   = mpatches.Patch(color="#E63946", label="Low coverage (<0.5x)")
    teal_patch  = mpatches.Patch(color="#2A9D8F", label="High coverage (>2.5x)")
    gray_patch  = mpatches.Patch(color="#aaaaaa", label="Normal coverage")
    ax.legend(handles=[red_patch, teal_patch, gray_patch],
              loc="upper right", fontsize=8)

    # ── Panel 2: Heatmap of top altered genes ──────────────
    ax = ax2

    # Select top 30 most interesting genes
    hits = df[df["any_alteration"]].copy()
    # Score: cancer_relevant > somatic_sv > cov_change
    hits["priority"] = (
        hits["cancer_relevant"].astype(int) * 4 +
        hits["somatic_sv"].astype(int) * 2 +
        hits["somatic_cov_change"].astype(int) * 1 +
        (hits["n_sv_pdx"] > 0).astype(int)
    )
    top30 = hits.nlargest(30, "priority")

    if len(top30) > 0:
        heatmap_data = top30.set_index("gene")[[
            "somatic_sv",
            "somatic_cov_change",
            "n_sv_pdx",
            "pdx_cov_ratio",
            "cancer_relevant",
        ]].copy()
        heatmap_data.columns = [
            "Somatic SV", "Somatic Cov\nChange",
            "# PDX SVs", "PDX Cov\nRatio", "Cancer\nRelevant"
        ]

        # Normalize for display
        display = heatmap_data.copy().astype(float)
        display["# PDX SVs"] = display["# PDX SVs"].clip(0, 5) / 5
        display["PDX Cov\nRatio"] = (display["PDX Cov\nRatio"] - 1).clip(-1, 2) / 2

        sns.heatmap(display, ax=ax, cmap="RdYlGn_r",
                    center=0, vmin=-0.5, vmax=1,
                    linewidths=0.5, linecolor="white",
                    cbar_kws={"shrink": 0.6, "label": "Scaled value"},
                    yticklabels=True)
        ax.set_title("Top Altered Dark Region Genes", fontsize=10, loc="left")
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", labelsize=8)
    else:
        ax.text(0.5, 0.5, "No altered genes detected\n(check thresholds)",
                ha="center", va="center", transform=ax.transAxes, fontsize=11)
        ax.set_title("Altered Dark Region Genes", fontsize=10)

    # ── Panel 3: Normal vs PDX coverage scatter ─────────────
    ax = ax3
    ax.set_facecolor("#f8f9fa")
    ax.spines[["top","right"]].set_visible(False)

    has_cov = df[df["normal_cov_ratio"].notna() & df["pdx_cov_ratio"].notna()]

    # Background: all genes
    ax.scatter(has_cov["normal_cov_ratio"], has_cov["pdx_cov_ratio"],
               c="#cccccc", s=8, alpha=0.4, label="All dark genes", zorder=1)

    # Highlight cancer-relevant
    cancer = has_cov[has_cov["cancer_relevant"]]
    ax.scatter(cancer["normal_cov_ratio"], cancer["pdx_cov_ratio"],
               c="#E63946", s=30, alpha=0.8, label="Cancer-relevant", zorder=3, marker="*")

    # Diagonal = no change
    lim = max(has_cov["normal_cov_ratio"].quantile(0.98),
              has_cov["pdx_cov_ratio"].quantile(0.98), 3.0)
    ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, alpha=0.5)
    ax.axhline(0.5, color="#E63946", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.axhline(2.5, color="#2A9D8F", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.axvline(0.5, color="#E63946", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.axvline(2.5, color="#2A9D8F", linewidth=0.5, linestyle=":", alpha=0.6)

    # Annotate genes that moved quadrants (somatic cov change)
    somatic_cov = has_cov[has_cov["somatic_cov_change"] & has_cov["cancer_relevant"]]
    for _, row in somatic_cov.iterrows():
        ax.annotate(row["gene"],
                    xy=(row["normal_cov_ratio"], row["pdx_cov_ratio"]),
                    fontsize=7, color="#9B2226", fontweight="bold",
                    xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("Normal Coverage Ratio", fontsize=10)
    ax.set_ylabel("PDX Coverage Ratio", fontsize=10)
    ax.set_title("Normal vs PDX Coverage\n(dark region genes)",
                 fontsize=10, loc="left")
    ax.set_xlim(-0.1, lim)
    ax.set_ylim(-0.1, lim)
    ax.legend(fontsize=8, loc="upper left")

    # ── Summary stats in figure ─────────────────────────────
    n_total   = len(df)
    n_altered = df["any_alteration"].sum()
    n_cancer  = (df["cancer_relevant"] & df["any_alteration"]).sum()
    n_somatic = df["somatic_sv"].sum()
    fig.text(0.01, 0.01,
             f"Dark region genes: {n_total}  |  "
             f"With any alteration: {n_altered}  |  "
             f"Cancer-relevant + altered: {n_cancer}  |  "
             f"Somatic SVs in dark regions: {n_somatic}",
             fontsize=8, color="#555", ha="left")

    plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--out", default="dark_region_results.png")
    args = parser.parse_args()
    df = load(args.report)
    plot_results(df, args.out)


if __name__ == "__main__":
    main()
