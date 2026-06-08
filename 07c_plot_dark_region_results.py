"""
07_plot_dark_region_results.py
------------------------------
Visualizes novel T2T gene alteration findings — all genes shown equally,
no pre-assumed cancer relevance filter.

  Panel 1 — Genome-wide coverage ratio per gene, all labelled if outlier
  Panel 2 — Heatmap of all altered genes (top 30 by evidence strength)
  Panel 3 — PDX coverage ratio distribution (all genes), or normal vs PDX scatter

Usage:
    python 07_plot_dark_region_results.py \
        --report dark_genes_report.tsv \
        --out    novel_t2t_results.png

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


def load(report_path):
    df = pd.read_csv(report_path, sep="\t")
    for col in ["normal_cov_ratio", "pdx_cov_ratio", "normal_depth", "pdx_depth"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df.attrs["has_normal"] = "normal_cov_ratio" in df.columns
    return df


def plot_results(df, out):
    has_normal = df.attrs.get("has_normal", False)

    fig = plt.figure(figsize=(18, 16))
    title_suffix = "PDX vs Normal" if has_normal else "PDX only"
    fig.suptitle(f"Novel T2T-CHM13 Gene Alterations — {title_suffix}",
                 fontsize=14, fontweight="bold", y=0.99)

    gs = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, :])   # full-width top panel
    ax2 = fig.add_subplot(gs[1, 0])   # heatmap
    ax3 = fig.add_subplot(gs[1, 1])   # distribution / scatter

    # ── Panel 1: Per-chromosome coverage ratio, ALL genes ──
    ax = ax1
    ax.set_facecolor("#f8f9fa")
    ax.spines[["top", "right"]].set_visible(False)

    chrom_order_present = [c for c in CHROM_ORDER if c in df["chrom"].values]
    x_pos = []
    x_labels = []
    cursor = 0

    for chrom in chrom_order_present:
        sub = df[df["chrom"] == chrom].copy()
        if len(sub) == 0:
            cursor += 0.5
            continue

        sub = sub.reset_index(drop=True)
        sub["x"] = cursor + np.linspace(0, 0.8, len(sub))
        x_pos.append(cursor + 0.4)
        x_labels.append(chrom.replace("chr", ""))

        for _, row in sub.iterrows():
            ratio = row["pdx_cov_ratio"]
            if pd.isna(ratio):
                continue

            is_outlier = ratio < 0.5 or ratio > 2.5
            has_sv     = row.get("n_sv_pdx", 0) > 0

            # Color by coverage status
            color = "#E63946" if ratio < 0.5 else \
                    "#2A9D8F" if ratio > 2.5 else \
                    "#aaaaaa"

            # Size/marker: larger if has SV, even larger if both SV + cov outlier
            if has_sv and is_outlier:
                size = 60; marker = "D"; alpha = 0.95   # diamond = both
            elif has_sv:
                size = 35; marker = "^"; alpha = 0.85   # triangle = SV only
            elif is_outlier:
                size = 30; marker = "o"; alpha = 0.85   # circle = cov only
            else:
                size = 10; marker = "o"; alpha = 0.4    # small = normal

            ax.scatter(row["x"], ratio, c=color, s=size,
                       alpha=alpha, marker=marker,
                       zorder=3 if (is_outlier or has_sv) else 2)

            # Label ALL outliers (coverage or SV)
            if is_outlier or has_sv:
                ax.annotate(
                    row["gene"],
                    xy=(row["x"], ratio),
                    xytext=(0, 8), textcoords="offset points",
                    fontsize=6, ha="center", color=color,
                    fontweight="bold",
                    rotation=45,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=color))

        cursor += 1.5

    ax.axhline(1.0, color="#333",    linewidth=0.8, linestyle="--", alpha=0.5)
    ax.axhline(0.5, color="#E63946", linewidth=0.6, linestyle=":",  alpha=0.7)
    ax.axhline(2.5, color="#2A9D8F", linewidth=0.6, linestyle=":",  alpha=0.7)
    ax.text(cursor + 0.1, 0.52, "del (<0.5x)", fontsize=7, color="#E63946")
    ax.text(cursor + 0.1, 2.52, "dup (>2.5x)", fontsize=7, color="#2A9D8F")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.set_ylabel("PDX Coverage Ratio\n(vs genome median)", fontsize=10)
    ax.set_title(
        "Novel T2T Gene Coverage per Chromosome  "
        "(◆ = SV + cov outlier,  ▲ = SV only,  ● = cov outlier)",
        fontsize=9, loc="left")
    ax.set_ylim(-0.1, max(4.0, df["pdx_cov_ratio"].dropna().quantile(0.99) + 0.5))

    red_patch  = mpatches.Patch(color="#E63946", label="Low coverage (<0.5x)")
    teal_patch = mpatches.Patch(color="#2A9D8F", label="High coverage (>2.5x)")
    gray_patch = mpatches.Patch(color="#aaaaaa", label="Normal coverage")
    ax.legend(handles=[red_patch, teal_patch, gray_patch],
              loc="upper right", fontsize=8)

    # ── Panel 2: Heatmap — all altered genes, top 30 by evidence ──
    ax = ax2

    hits = df[df["any_alteration"]].copy()

    # Priority score — evidence-based only, no cancer assumption
    hits["priority"] = (hits["n_sv_pdx"] > 0).astype(int) * 2
    hits["priority"] += hits["pdx_cov_ratio"].apply(
        lambda r: 2 if (not pd.isna(r) and (r < 0.5 or r > 2.5)) else 0)
    if has_normal and "somatic_sv" in hits.columns:
        hits["priority"] += hits["somatic_sv"].astype(int) * 3
    if has_normal and "somatic_cov_change" in hits.columns:
        hits["priority"] += hits["somatic_cov_change"].astype(int) * 2

    top30 = hits.nlargest(30, "priority")

    if len(top30) > 0:
        hmap_cols   = []
        hmap_labels = []
        if has_normal and "somatic_sv" in top30.columns:
            hmap_cols.append("somatic_sv")
            hmap_labels.append("Somatic SV")
        if has_normal and "somatic_cov_change" in top30.columns:
            hmap_cols.append("somatic_cov_change")
            hmap_labels.append("Somatic Cov\nChange")
        hmap_cols   += ["n_sv_pdx", "pdx_cov_ratio"]
        hmap_labels += ["# PDX SVs", "PDX Cov\nRatio"]

        heatmap_data = top30.set_index("gene")[hmap_cols].copy()
        heatmap_data.columns = hmap_labels

        display = heatmap_data.copy().astype(float)
        display["# PDX SVs"]     = display["# PDX SVs"].clip(0, 5) / 5
        display["PDX Cov\nRatio"] = (display["PDX Cov\nRatio"] - 1).clip(-1, 2) / 2

        sns.heatmap(display, ax=ax, cmap="RdYlGn_r",
                    center=0, vmin=-0.5, vmax=1,
                    linewidths=0.5, linecolor="white",
                    cbar_kws={"shrink": 0.6, "label": "Scaled value"},
                    yticklabels=True)
        ax.set_title("Top 30 Altered Novel T2T Genes\n(ranked by evidence strength)",
                     fontsize=10, loc="left")
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", labelsize=8)
    else:
        ax.text(0.5, 0.5, "No altered genes detected\n(check thresholds)",
                ha="center", va="center", transform=ax.transAxes, fontsize=11)
        ax.set_title("Altered Novel T2T Genes", fontsize=10)

    # ── Panel 3: Coverage distribution or normal vs PDX scatter ──
    ax = ax3
    ax.set_facecolor("#f8f9fa")
    ax.spines[["top", "right"]].set_visible(False)

    if has_normal and "normal_cov_ratio" in df.columns:
        # Normal vs PDX scatter — all genes shown
        has_cov = df[df["normal_cov_ratio"].notna() & df["pdx_cov_ratio"].notna()]
        ax.scatter(has_cov["normal_cov_ratio"], has_cov["pdx_cov_ratio"],
                   c="#aaaaaa", s=10, alpha=0.5, label="All novel genes", zorder=1)
        # Highlight any outlier (moved quadrant)
        outliers = has_cov[has_cov["any_alteration"]]
        ax.scatter(outliers["normal_cov_ratio"], outliers["pdx_cov_ratio"],
                   c="#E63946", s=25, alpha=0.85, label="Altered", zorder=3)
        for _, row in outliers.iterrows():
            ax.annotate(row["gene"],
                        xy=(row["normal_cov_ratio"], row["pdx_cov_ratio"]),
                        fontsize=6.5, color="#9B2226", fontweight="bold",
                        xytext=(4, 4), textcoords="offset points")
        lim = max(has_cov["normal_cov_ratio"].quantile(0.98),
                  has_cov["pdx_cov_ratio"].quantile(0.98), 3.0)
        ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, alpha=0.5)
        ax.axhline(0.5, color="#E63946", linewidth=0.5, linestyle=":", alpha=0.6)
        ax.axhline(2.5, color="#2A9D8F", linewidth=0.5, linestyle=":", alpha=0.6)
        ax.axvline(0.5, color="#E63946", linewidth=0.5, linestyle=":", alpha=0.6)
        ax.axvline(2.5, color="#2A9D8F", linewidth=0.5, linestyle=":", alpha=0.6)
        ax.set_xlabel("Normal Coverage Ratio", fontsize=10)
        ax.set_ylabel("PDX Coverage Ratio", fontsize=10)
        ax.set_title("Normal vs PDX Coverage\n(all novel T2T genes)",
                     fontsize=10, loc="left")
        ax.set_xlim(-0.1, lim); ax.set_ylim(-0.1, lim)
        ax.legend(fontsize=8, loc="upper left")

    else:
        # PDX-only: ranked coverage ratio — all genes, label all outliers
        has_cov    = df[df["pdx_cov_ratio"].notna()].copy()
        sorted_df  = has_cov.sort_values("pdx_cov_ratio").reset_index(drop=True)
        ratios     = sorted_df["pdx_cov_ratio"].values

        colors = ["#E63946" if r < 0.5 else
                  "#2A9D8F" if r > 2.5 else
                  "#aaaaaa" for r in ratios]

        ax.scatter(range(len(sorted_df)), ratios,
                   c=colors, s=10, alpha=0.6, zorder=2)

        # Label every outlier gene
        for i, row in sorted_df[
                (sorted_df["pdx_cov_ratio"] < 0.5) |
                (sorted_df["pdx_cov_ratio"] > 2.5)].iterrows():
            color = "#E63946" if row["pdx_cov_ratio"] < 0.5 else "#2A9D8F"
            ax.annotate(row["gene"],
                        xy=(i, row["pdx_cov_ratio"]),
                        xytext=(0, 6), textcoords="offset points",
                        fontsize=6.5, ha="center",
                        color=color, fontweight="bold", rotation=45,
                        arrowprops=dict(arrowstyle="-", lw=0.3, color=color))

        ax.axhline(1.0, color="#333",    linewidth=0.8, linestyle="--", alpha=0.5)
        ax.axhline(0.5, color="#E63946", linewidth=0.6, linestyle=":",  alpha=0.7)
        ax.axhline(2.5, color="#2A9D8F", linewidth=0.6, linestyle=":",  alpha=0.7)
        ax.set_xlabel("Genes (ranked by coverage ratio)", fontsize=10)
        ax.set_ylabel("PDX Coverage Ratio", fontsize=10)
        ax.set_title("PDX Coverage Ratio — All Novel T2T Genes\n(all outliers labelled)",
                     fontsize=10, loc="left")

        red_p  = mpatches.Patch(color="#E63946", label=f"Low (<0.5x):  {(ratios<0.5).sum()} genes")
        teal_p = mpatches.Patch(color="#2A9D8F", label=f"High (>2.5x): {(ratios>2.5).sum()} genes")
        gray_p = mpatches.Patch(color="#aaaaaa", label=f"Normal:       {((ratios>=0.5)&(ratios<=2.5)).sum()} genes")
        ax.legend(handles=[red_p, teal_p, gray_p], fontsize=8, loc="upper left")

    # ── Footer summary ──────────────────────────────────────
    n_total   = len(df)
    n_altered = df["any_alteration"].sum()
    n_sv      = (df["n_sv_pdx"] > 0).sum()
    n_low_cov = (df["pdx_cov_ratio"] < 0.5).sum()
    n_hi_cov  = (df["pdx_cov_ratio"] > 2.5).sum()
    summary   = (f"Novel T2T genes: {n_total}  |  "
                 f"Any alteration: {n_altered}  |  "
                 f"With SVs: {n_sv}  |  "
                 f"Low coverage: {n_low_cov}  |  "
                 f"High coverage: {n_hi_cov}")
    if has_normal and "somatic_sv" in df.columns:
        summary += f"  |  Somatic SVs: {df['somatic_sv'].sum()}"
    fig.text(0.01, 0.005, summary, fontsize=8, color="#555", ha="left")

    plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True,
                        help="dark_genes_report.tsv from 06_intersect_annotate.py")
    parser.add_argument("--out", default="novel_t2t_results.png")
    args = parser.parse_args()
    df = load(args.report)
    plot_results(df, args.out)


if __name__ == "__main__":
    main()
