"""
02_extract_dark_genes.py
------------------------
Parses the CHM13 v2.0 GFF3 to extract protein-coding genes,
then intersects with segmental duplications + centromere satellites
to define "dark region genes" — genes that fall in regions
inaccessible or misassembled in hg38.

Also flags genes with known cancer relevance.

Run AFTER 01_download_chm13_annotation.sh

Usage:
    python 02_extract_dark_genes.py \
        --gff   chm13_annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz \
        --segdup chm13_annotation/chm13v2.0_segmental_duplications.bed \
        --censat chm13_annotation/chm13v2.0_censat_v2.1.bed \
        --out_bed dark_region_genes.bed \
        --out_all all_chm13_genes.bed

Dependencies:
    pip install --user pandas
    module load bedtools
"""

import argparse
import gzip
import subprocess
import os
import pandas as pd

# ── Known cancer-relevant genes in T2T-newly-resolved regions ──────────────
# Many of these were unmappable/collapsed in hg38
CANCER_GENES = {
    # Prostate cancer specific
    "KLK3","KLK2","KLK4","AR","FOXA1","ERG","PTEN","TP53","RB1",
    "BRCA2","CDK12","SPOP","MED12","CHD1","NCOR1","ZBTB16",
    # Dark region / segdup genes with cancer links
    "NOTCH2NL","GOLGA6L1","GOLGA6L2","GOLGA6L4","GOLGA8A","GOLGA8B",
    "NBPF10","NBPF14","NBPF15","NBPF19","NBPF20",  # neuroblastoma amplified
    "SPDYE1","SPDYE2","SPDYE3","SPDYE4","SPDYE5","SPDYE6",
    "SMN1","SMN2",          # chr5 dark region
    "STRC","OTOA",          # chr15/16 segdup
    "CFC1B",                # chr2 dark
    "RBMXL1","RBMXL2","RBMXL3",  # chr15 dark
    "ARHGAP11A","ARHGAP11B",      # chr15 segdup
    "SERF1A","SERF1B",
    "GTF2I","GTF2IRD1",    # chr7 segdup
    "HYDIN","HYDIN2",       # chr16 dark
    "CHRNA7","CHRFAM7A",    # chr15
    "CFHR1","CFHR3",        # chr1 segdup
    "C4A","C4B",            # chr6 MHC dark
    "GBA","GBA2",           # chr1
    "PMS2","PMS2CL",        # mismatch repair
}

def parse_gff3(gff_path):
    """Parse GFF3 (gzipped or not) → DataFrame of gene records."""
    genes = []
    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attrs = cols
            if feature != "gene":
                continue
            # Parse attributes
            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k.strip()] = v.strip()
            gene_name = attr_dict.get("gene", attr_dict.get("Name", attr_dict.get("ID", ".")))
            gene_type  = attr_dict.get("gene_biotype", attr_dict.get("biotype", "unknown"))
            genes.append({
                "chrom":  seqname,
                "start":  int(start) - 1,   # GFF3 is 1-based, convert to 0-based
                "end":    int(end),
                "name":   gene_name,
                "biotype": gene_type,
                "strand": strand,
            })
    return pd.DataFrame(genes)


def intersect_with_bed(genes_bed_path, region_bed_path, label):
    """Use bedtools intersect to find genes overlapping a region BED."""
    cmd = [
        "bedtools", "intersect",
        "-a", genes_bed_path,
        "-b", region_bed_path,
        "-wa",          # report original gene record
        "-u",           # only once even if multiple overlaps
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  bedtools error for {label}: {result.stderr[:200]}")
        return set()
    names = set()
    for line in result.stdout.strip().split("\n"):
        if line:
            cols = line.split("\t")
            if len(cols) >= 4:
                names.add(cols[3])
    print(f"  {label}: {len(names)} genes overlap")
    return names


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff",     required=True)
    parser.add_argument("--segdup",  required=True)
    parser.add_argument("--censat",  required=True)
    parser.add_argument("--out_bed", default="dark_region_genes.bed")
    parser.add_argument("--out_all", default="all_chm13_genes.bed")
    args = parser.parse_args()

    print("Parsing GFF3...")
    df = parse_gff3(args.gff)
    print(f"  Total gene records: {len(df)}")

    # Keep protein-coding + lncRNA (most relevant)
    coding = df[df["biotype"].isin(["protein_coding","lncRNA","miRNA","snRNA","snoRNA"])]
    print(f"  Protein-coding / ncRNA: {len(coding)}")

    # Save all genes BED (used later for SV intersection)
    all_bed = args.out_all
    coding[["chrom","start","end","name","biotype","strand"]].to_csv(
        all_bed, sep="\t", header=False, index=False)
    print(f"Saved all genes BED: {all_bed}")

    # Write temp BED for intersect
    tmp_bed = "/tmp/chm13_genes_tmp.bed"
    coding[["chrom","start","end","name","biotype","strand"]].to_csv(
        tmp_bed, sep="\t", header=False, index=False)

    # Intersect with dark regions
    print("\nIntersecting with dark region tracks...")
    segdup_genes = intersect_with_bed(tmp_bed, args.segdup, "Segmental duplications")
    censat_genes = intersect_with_bed(tmp_bed, args.censat, "Centromere/satellites")

    dark_names = segdup_genes | censat_genes
    print(f"\n  Total dark region genes: {len(dark_names)}")

    # Flag cancer relevance
    coding = coding.copy()
    coding["dark_region"] = coding["name"].isin(dark_names)
    coding["cancer_relevant"] = coding["name"].isin(CANCER_GENES)
    coding["dark_and_cancer"] = coding["dark_region"] & coding["cancer_relevant"]

    # Save dark region gene BED
    dark_df = coding[coding["dark_region"]].copy()
    dark_df["score"] = dark_df["cancer_relevant"].astype(int)  # 1=cancer, 0=not
    dark_df[["chrom","start","end","name","score","strand"]].to_csv(
        args.out_bed, sep="\t", header=False, index=False)
    print(f"Saved dark region gene BED: {args.out_bed}")

    # Summary
    print(f"\n{'='*50}")
    print(f"SUMMARY")
    print(f"{'='*50}")
    print(f"All CHM13 genes parsed:           {len(df):>6}")
    print(f"Protein-coding / ncRNA:           {len(coding):>6}")
    print(f"In dark regions (segdup/censat):  {len(dark_df):>6}")
    print(f"Dark + cancer-relevant:           {coding['dark_and_cancer'].sum():>6}")

    print(f"\nCancer-relevant dark region genes:")
    cancer_dark = coding[coding["dark_and_cancer"]].sort_values("chrom")
    print(cancer_dark[["chrom","start","end","name","biotype"]].to_string(index=False))

    # Save full annotated table
    annotated_out = args.out_bed.replace(".bed","_annotated.tsv")
    coding.to_csv(annotated_out, sep="\t", index=False)
    print(f"\nFull annotated table: {annotated_out}")


if __name__ == "__main__":
    main()
