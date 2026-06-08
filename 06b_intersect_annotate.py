"""
06_intersect_annotate.py
------------------------
Integrates three evidence layers for each dark region gene:
  1. SV overlap  — does an SV breakpoint fall in or near this gene?
  2. Coverage    — is depth significantly higher/lower than genome median?
  3. Somatic     — is the event PDX-specific (not in normal)?

Outputs:
  - dark_genes_report.tsv   — full annotated table
  - dark_genes_hits.tsv     — genes with any evidence of alteration
  - (used by 07_plot_dark_regions.py)

Usage:
    python 06_intersect_annotate.py \
        --dark_bed      dark_region_genes.bed \
        --normal_vcf    vcf_chm13/normal_chm13_germline.vcf \
        --pdx_vcf       vcf_chm13/PDX_chm13_germline.vcf \
        --normal_cov    mosdepth_out/normal_dark.regions.bed \
        --pdx_cov       mosdepth_out/PDX_dark.regions.bed \
        --normal_genome mosdepth_out/normal_genome.cov \
        --pdx_genome    mosdepth_out/PDX_genome.cov \

    --normal_cov / --pdx_cov:
        Output of coverage_samtools_depth.py  (5-col BED: chrom start end name depth)
        Also accepts mosdepth .regions.bed.gz if available

    --normal_genome / --pdx_genome:
        Output of: samtools coverage --min-MQ 1 your.bam > sample_genome.cov
        Also accepts mosdepth genome windows .bed.gz if available\
        --out_report    dark_genes_report.tsv \
        --out_hits      dark_genes_hits.tsv

Dependencies:
    pip install --user pandas numpy scipy
    module load bedtools
"""

import argparse
import gzip
import subprocess
import re
import numpy as np
import pandas as pd
from scipy import stats

# ── Cancer gene list (same as step 2) ──────────────────────
CANCER_GENES = {
    "KLK3","KLK2","KLK4","AR","FOXA1","ERG","PTEN","TP53","RB1",
    "BRCA2","CDK12","SPOP","MED12","CHD1","NCOR1","ZBTB16",
    "NOTCH2NL","GOLGA6L1","GOLGA6L2","GOLGA6L4","GOLGA8A","GOLGA8B",
    "NBPF10","NBPF14","NBPF15","NBPF19","NBPF20",
    "SPDYE1","SPDYE2","SPDYE3","SPDYE4","SPDYE5","SPDYE6",
    "SMN1","SMN2","STRC","OTOA","CFC1B",
    "RBMXL1","RBMXL2","RBMXL3",
    "ARHGAP11A","ARHGAP11B","SERF1A","SERF1B",
    "GTF2I","GTF2IRD1","HYDIN","HYDIN2","CHRNA7","CHRFAM7A",
    "CFHR1","CFHR3","C4A","C4B","GBA","GBA2","PMS2","PMS2CL",
}

# ── Thresholds ──────────────────────────────────────────────
COV_LOW_THRESH  = 0.5   # depth < 50% of genome median → potential deletion
COV_HIGH_THRESH = 2.5   # depth > 2.5x genome median → potential amplification
FLANK_BP        = 500   # bp flanking gene to catch nearby breakpoints


def load_bed(path):
    """Load BED file (gzipped or not) into DataFrame."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 4:
                rows.append(cols)
    if not rows:
        return pd.DataFrame(columns=["chrom","start","end","name"])
    df = pd.DataFrame(rows)
    df.columns = range(len(df.columns))
    return df


def parse_vcf_to_bed(vcf_path, flank=FLANK_BP):
    """Parse VCF → BED of SV breakpoints (PASS only), with flanking."""
    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, sv_id, ref, alt, qual, filt, info = cols[:8]
            if filt != "PASS":
                continue
            svtype = re.search(r'SVTYPE=([^;]+)', info)
            svtype = svtype.group(1) if svtype else "UNK"
            svlen  = re.search(r'SVLEN=(-?\d+)', info)
            svlen  = abs(int(svlen.group(1))) if svlen else 0
            end_m  = re.search(r'END=(\d+)', info)
            end    = int(end_m.group(1)) if end_m else int(pos)
            support = re.search(r'SUPPORT=(\d+)', info)
            support = int(support.group(1)) if support else 0
            vaf    = re.search(r'VAF=([\d.]+)', info)
            vaf    = float(vaf.group(1)) if vaf else 0.0
            records.append({
                "chrom": chrom,
                "start": max(0, int(pos) - flank),
                "end":   end + flank,
                "sv_id": sv_id,
                "svtype": svtype,
                "svlen": svlen,
                "support": support,
                "vaf": vaf,
            })
    return pd.DataFrame(records)


def bedtools_intersect(a_df, b_df, a_name="a", b_name="b"):
    """
    Run bedtools intersect in memory via temp files.
    Returns set of (chrom, start, end, name) tuples from A that overlap B.
    """
    import tempfile, os
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as fa:
        for _, r in a_df.iterrows():
            fa.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r.get('name','.')}\n")
        fa_path = fa.name
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as fb:
        for _, r in b_df.iterrows():
            fb.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r.get('sv_id', r.get('name','.'))}\n")
        fb_path = fb.name

    result = subprocess.run(
        ["bedtools", "intersect", "-a", fa_path, "-b", fb_path, "-wa", "-wb"],
        capture_output=True, text=True
    )
    os.unlink(fa_path); os.unlink(fb_path)

    hits = {}  # gene_name -> list of SV records
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        gene_name = cols[3]
        sv_name   = cols[7] if len(cols) > 7 else "."
        if gene_name not in hits:
            hits[gene_name] = []
        hits[gene_name].append(sv_name)
    return hits


def genome_median_depth(genome_cov_path):
    """
    Compute genome-wide median depth.
    Accepts two formats:
      1. samtools coverage output  (tab-separated, header line starts with #rname)
         columns: #rname startpos endpos numreads covbases coverage meandepth ...
      2. mosdepth genome windows   (.bed or .bed.gz)
         columns: chrom start end depth
    Auto-detects format from file header.
    """
    opener = gzip.open if genome_cov_path.endswith(".gz") else open
    depths = []
    lengths = []
    is_samtools_cov = False

    with opener(genome_cov_path, "rt") as f:
        for line in f:
            # samtools coverage header
            if line.startswith("#rname") or line.startswith("#"):
                if "meandepth" in line:
                    is_samtools_cov = True
                continue
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue

            if is_samtools_cov:
                # samtools coverage: col6=meandepth, col2-col1=length
                try:
                    length    = int(cols[2]) - int(cols[1])
                    meandepth = float(cols[6])
                    # Skip unplaced/small contigs
                    chrom = cols[0]
                    if length < 1_000_000:
                        continue
                    depths.append(meandepth)
                    lengths.append(length)
                except (ValueError, IndexError):
                    continue
            else:
                # mosdepth windows BED: col3=depth
                try:
                    depths.append(float(cols[3]))
                    lengths.append(1)
                except (ValueError, IndexError):
                    continue

    if not depths:
        print(f"  WARNING: could not parse genome coverage from {genome_cov_path}, using 30x default")
        return 30.0

    if is_samtools_cov:
        # Weighted median by chromosome length
        total = sum(lengths)
        pairs = sorted(zip(depths, lengths), key=lambda x: x[0])
        cumlen = 0
        for d, l in pairs:
            cumlen += l
            if cumlen >= total / 2:
                return float(d)
        return float(pairs[-1][0])
    else:
        return float(np.median(depths))


def load_region_coverage(cov_path):
    """
    Load per-gene coverage into dict of gene_name -> depth.
    Accepts:
      1. coverage_samtools_depth.py output  (plain .bed, col4=name, col5=depth)
      2. mosdepth .regions.bed.gz           (col4=name, col5=depth)
    Both have identical format — 5 columns: chrom start end name depth.
    """
    opener = gzip.open if cov_path.endswith(".gz") else open
    cov = {}
    with opener(cov_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue
            name  = cols[3]
            try:
                depth = float(cols[4])
            except ValueError:
                continue
            cov[name] = depth
    return cov


def main():
    parser = argparse.ArgumentParser(
        description="Intersect SVs + coverage with dark region genes",
        formatter_class=argparse.RawTextHelpFormatter)

    # Always required
    parser.add_argument("--dark_bed",   required=True,
                        help="Dark region gene BED file")
    parser.add_argument("--pdx_vcf",    required=True,
                        help="PDX Sniffles2 VCF (PASS SVs)")
    parser.add_argument("--pdx_cov",    required=True,
                        help="PDX per-gene coverage BED (from coverage_samtools_depth.py)")
    parser.add_argument("--pdx_genome", required=True,
                        help="PDX genome-wide coverage (samtools coverage output)")

    # Optional normal tissue inputs — if omitted, somatic comparisons are skipped
    parser.add_argument("--normal_vcf",    default=None,
                        help="[optional] Normal tissue VCF — enables somatic SV detection")
    parser.add_argument("--normal_cov",    default=None,
                        help="[optional] Normal per-gene coverage BED — enables somatic coverage change")
    parser.add_argument("--normal_genome", default=None,
                        help="[optional] Normal genome-wide coverage (samtools coverage output)")

    parser.add_argument("--out_report", default="dark_genes_report.tsv")
    parser.add_argument("--out_hits",   default="dark_genes_hits.tsv")
    args = parser.parse_args()

    has_normal = args.normal_vcf is not None
    if not has_normal:
        print("NOTE: No normal tissue inputs provided — somatic comparison disabled.")
        print("      Reporting PDX-only SV and coverage alterations.\n")

    # ── Load dark gene BED ─────────────────────────────────
    print("Loading dark region genes...")
    dark_df = load_bed(args.dark_bed)
    dark_df.columns = ["chrom","start","end","name","score","strand"] + \
                      list(range(len(dark_df.columns)-6))
    dark_df["start"] = dark_df["start"].astype(int)
    dark_df["end"]   = dark_df["end"].astype(int)
    dark_df["cancer_relevant"] = dark_df["name"].isin(CANCER_GENES)
    print(f"  {len(dark_df)} dark region genes loaded")

    # ── Parse SVs ─────────────────────────────────────────
    print("Parsing VCFs...")
    pdx_svs = parse_vcf_to_bed(args.pdx_vcf)
    print(f"  PDX SVs (PASS): {len(pdx_svs)}")

    normal_svs  = pd.DataFrame()
    normal_hits = {}
    if has_normal:
        normal_svs = parse_vcf_to_bed(args.normal_vcf)
        print(f"  Normal SVs (PASS): {len(normal_svs)}")

    # ── Intersect SVs with dark genes ─────────────────────
    print("Intersecting SVs with dark region genes...")
    pdx_hits = bedtools_intersect(dark_df, pdx_svs)
    print(f"  Genes hit in PDX: {len(pdx_hits)}")

    if has_normal:
        normal_hits = bedtools_intersect(dark_df, normal_svs)
        print(f"  Genes hit in normal: {len(normal_hits)}")

    # ── Load coverage ──────────────────────────────────────
    print("Loading coverage data...")
    pdx_med = genome_median_depth(args.pdx_genome)
    pdx_cov = load_region_coverage(args.pdx_cov)
    print(f"  PDX genome median depth: {pdx_med:.1f}x")

    normal_med = np.nan
    normal_cov = {}
    if has_normal and args.normal_cov and args.normal_genome:
        normal_med = genome_median_depth(args.normal_genome)
        normal_cov = load_region_coverage(args.normal_cov)
        print(f"  Normal genome median depth: {normal_med:.1f}x")

    # ── Build report ───────────────────────────────────────
    print("Building report...")

    def cov_class(ratio):
        if np.isnan(ratio):          return "no_data"
        if ratio < COV_LOW_THRESH:   return "low_coverage"
        if ratio > COV_HIGH_THRESH:  return "high_coverage"
        return "normal"

    rows = []
    for _, gene in dark_df.iterrows():
        gname = gene["name"]

        # SV evidence
        n_sv_pdx   = len(pdx_hits.get(gname, []))
        sv_ids_pdx = ";".join(pdx_hits.get(gname, [])) or "."
        n_sv_normal = len(normal_hits.get(gname, [])) if has_normal else None
        # Somatic SV = in PDX but not in normal (only if normal provided)
        somatic_sv = (n_sv_pdx > 0 and n_sv_normal == 0) if has_normal else None

        # PDX coverage
        pdx_depth = pdx_cov.get(gname, np.nan)
        pdx_ratio = pdx_depth / pdx_med if not np.isnan(pdx_depth) and pdx_med > 0 else np.nan
        pdx_cov_class = cov_class(pdx_ratio)

        # Normal coverage (optional)
        norm_depth = normal_cov.get(gname, np.nan) if has_normal else np.nan
        norm_ratio = norm_depth / normal_med if (
            has_normal and not np.isnan(norm_depth) and not np.isnan(normal_med)
            and normal_med > 0) else np.nan
        norm_cov_class = cov_class(norm_ratio) if has_normal else "not_run"

        # Somatic coverage change = normal OK, PDX abnormal
        somatic_cov_change = (
            norm_cov_class == "normal" and
            pdx_cov_class in ("low_coverage", "high_coverage")
        ) if has_normal else None

        # any_hit: when no normal, flag based on PDX coverage alone
        if has_normal:
            any_hit = (n_sv_pdx > 0 or somatic_sv or
                       somatic_cov_change or pdx_cov_class != "normal")
        else:
            any_hit = (n_sv_pdx > 0 or pdx_cov_class != "normal")

        row = {
            "chrom":           gene["chrom"],
            "start":           gene["start"],
            "end":             gene["end"],
            "gene":            gname,
            "strand":          gene["strand"],
            "cancer_relevant": gene["cancer_relevant"],
            # PDX SV
            "n_sv_pdx":        n_sv_pdx,
            "pdx_sv_ids":      sv_ids_pdx,
            # PDX coverage
            "pdx_depth":       round(pdx_depth, 2) if not np.isnan(pdx_depth) else "NA",
            "pdx_cov_ratio":   round(pdx_ratio,  3) if not np.isnan(pdx_ratio)  else "NA",
            "pdx_cov_class":   pdx_cov_class,
            # Summary
            "any_alteration":  any_hit,
        }

        # Add normal columns only if normal provided
        if has_normal:
            row.update({
                "n_sv_normal":        n_sv_normal,
                "somatic_sv":         somatic_sv,
                "normal_depth":       round(norm_depth, 2) if not np.isnan(norm_depth) else "NA",
                "normal_cov_ratio":   round(norm_ratio,  3) if not np.isnan(norm_ratio)  else "NA",
                "normal_cov_class":   norm_cov_class,
                "somatic_cov_change": somatic_cov_change,
            })

        rows.append(row)

    report = pd.DataFrame(rows)

    # Sort: cancer-relevant first, then by evidence strength
    sort_cols = ["cancer_relevant", "any_alteration", "n_sv_pdx", "pdx_cov_ratio"]
    sort_asc  = [False, False, False, True]
    if has_normal:
        sort_cols = ["cancer_relevant", "somatic_sv",
                     "somatic_cov_change", "any_alteration", "n_sv_pdx"]
        sort_asc  = [False, False, False, False, False]
    report = report.sort_values(sort_cols, ascending=sort_asc)

    report.to_csv(args.out_report, sep="\t", index=False)
    print(f"\nFull report: {args.out_report}  ({len(report)} genes)")

    # ── Hits summary ───────────────────────────────────────
    hits = report[report["any_alteration"]].copy()
    hits.to_csv(args.out_hits, sep="\t", index=False)
    print(f"Genes with alterations: {args.out_hits}  ({len(hits)} genes)")

    # ── Print top findings ─────────────────────────────────
    print(f"\n{'='*60}")
    print("TOP FINDINGS — Cancer-relevant dark region genes altered in PDX")
    print(f"{'='*60}")
    top = hits[hits["cancer_relevant"]].head(20)
    if len(top) == 0:
        print("  (none — check all dark region hits in out_hits file)")
        top = hits.head(20)

    # Show columns that exist (normal cols may be absent)
    show_cols = ["gene","chrom","n_sv_pdx","pdx_cov_class","pdx_cov_ratio"]
    if has_normal:
        show_cols += ["somatic_sv","somatic_cov_change"]
    show_cols = [c for c in show_cols if c in top.columns]
    print(top[show_cols].to_string(index=False))


if __name__ == "__main__":
    main()
