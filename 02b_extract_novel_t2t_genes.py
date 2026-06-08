"""
02b_extract_novel_t2t_genes.py
------------------------------
Extracts genes present in T2T-CHM13 v5.2 but ABSENT from hg38 —
the ~99-140 newly annotated protein-coding genes.

Uses the JHU RefSeqv110 + Liftoff v5.2 GFF3 which already incorporates
CAT/Liftoff annotations. Novel genes are identified by:

  1. extra_copy_number > 0  → extra copy found only in CHM13 (not in hg38)
  2. gene_name starts with LOC  → NCBI ID for novel/unnamed gene
  3. Both together              → highest confidence novel call

NOTE: extra_copy_number=0 with a LOC name = primary copy of a gene that
exists in CHM13 but has no standard name in hg38. These are also novel.

Usage:
    python 02b_extract_novel_t2t_genes.py \
        --gff     chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz \
        --out     novel_t2t_genes.bed \
        --out_tsv novel_t2t_genes_annotated.tsv

    # Stricter: only protein-coding
    python 02b_extract_novel_t2t_genes.py \
        --gff     chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz \
        --out     novel_t2t_genes.bed \
        --coding_only

Dependencies:
    pip install --user pandas
"""

import argparse
import gzip
import re
import pandas as pd

# Cancer-relevant genes among T2T novel/extra-copy set
CANCER_RELEVANT = {
    "NOTCH2NL","NOTCH2NLA","NOTCH2NLB","NOTCH2NLC","NOTCH2NLD",
    "NBPF10","NBPF14","NBPF15","NBPF19","NBPF20","NBPF25P",
    "GOLGA6L1","GOLGA6L2","GOLGA6L4","GOLGA6L9","GOLGA6L10",
    "GOLGA8A","GOLGA8B","GOLGA8C","GOLGA8D","GOLGA8E",
    "WASH1","WASHC1","GPRIN2","LOC124900631",
    "SPDYE1","SPDYE2","SPDYE3","SPDYE4","SPDYE5","SPDYE6",
    "SMN1","SMN2","STRC","OTOA",
    "ARHGAP11B","SERF1A","SERF1B",
    "GTF2I","GTF2IRD1","HYDIN2","CHRFAM7A",
    "PMS2CL","GBA2","C4B_2",
}


def parse_gff3(gff_path, coding_only=False):
    """
    Parse CHM13 v5.2 GFF3 and extract novel genes using:
      - extra_copy_number > 0  (CHM13-specific extra copy)
      - gene_name = LOC*       (novel gene with no hg38 name)
    """
    genes = []
    opener = gzip.open if gff_path.endswith(".gz") else open

    with opener(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue

            seqname = cols[0]
            start   = cols[3]
            end     = cols[4]
            strand  = cols[6]
            attrs   = cols[8]

            # Parse attributes
            attr = {}
            for item in attrs.split(";"):
                item = item.strip()
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr[k.strip()] = v.strip()

            gene_name  = attr.get("gene_name", attr.get("gene", attr.get("Name", ".")))
            gene_type  = attr.get("gene_biotype", "unknown")
            extra_copy = int(attr.get("extra_copy_number", 0))
            coverage   = float(attr.get("coverage", 1.0))
            seq_id     = float(attr.get("sequence_ID", 1.0))

            # Skip pseudogenes
            if "pseudogene" in gene_type.lower():
                continue

            if coding_only and gene_type != "protein_coding":
                continue

            # Novel gene detection
            is_loc        = bool(re.match(r"LOC\d+", gene_name))
            is_extra_copy = extra_copy > 0
            is_novel      = is_loc or is_extra_copy

            if not is_novel:
                continue

            # Classify novelty type
            if is_loc and not is_extra_copy:
                novel_type = "novel_gene"        # entirely new, no hg38 equivalent
            elif is_extra_copy and not is_loc:
                novel_type = "extra_copy"         # known gene, CHM13-only extra copy
            else:
                novel_type = "novel_extra_copy"   # novel gene AND extra copy

            genes.append({
                "chrom":           seqname,
                "start":           int(start) - 1,   # GFF3 1-based -> 0-based BED
                "end":             int(end),
                "name":            gene_name,
                "gene_type":       gene_type,
                "strand":          strand,
                "extra_copy_num":  extra_copy,
                "coverage":        coverage,
                "sequence_id":     seq_id,
                "novel_type":      novel_type,
                "cancer_relevant": gene_name in CANCER_RELEVANT,
            })

    return pd.DataFrame(genes)


def main():
    parser = argparse.ArgumentParser(
        description="Extract T2T-CHM13 novel genes from v5.2 GFF3")
    parser.add_argument("--gff",         required=True,
                        help="chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz")
    parser.add_argument("--out",         default="novel_t2t_genes.bed",
                        help="Output BED for 06_intersect_annotate.py")
    parser.add_argument("--out_tsv",     default="novel_t2t_genes_annotated.tsv",
                        help="Full annotated TSV with metadata")
    parser.add_argument("--coding_only", action="store_true",
                        help="Keep only protein_coding genes (stricter)")
    parser.add_argument("--min_seq_id",  type=float, default=0.0,
                        help="Min sequence identity (0-1). 0=keep all, "
                             "0.8=filter low-quality copies. Default: 0")
    args = parser.parse_args()

    print(f"Parsing: {args.gff}")
    df = parse_gff3(args.gff, coding_only=args.coding_only)

    if args.min_seq_id > 0:
        before = len(df)
        df = df[df["sequence_id"] >= args.min_seq_id]
        print(f"  Filtered by sequence_id >= {args.min_seq_id}: {before} -> {len(df)}")

    if len(df) == 0:
        print("\nWARNING: No novel genes found.")
        print("Check GFF3 attributes with:")
        print("  zcat your.gff3.gz | awk '$3==\"gene\"' | head -3 | cut -f9 | tr ';' '\\n'")
        return

    # ── Summary ─────────────────────────────────────────────
    print(f"\n{'='*55}")
    print("NOVEL T2T GENE SUMMARY")
    print(f"{'='*55}")
    print(f"Total novel genes:              {len(df)}")
    print(f"  Novel (LOC name, new gene):   {(df['novel_type']=='novel_gene').sum()}")
    print(f"  Extra copy (known gene name): {(df['novel_type']=='extra_copy').sum()}")
    print(f"  Novel + extra copy:           {(df['novel_type']=='novel_extra_copy').sum()}")
    print(f"\nBiotype breakdown:")
    print(df["gene_type"].value_counts().to_string())
    print(f"\nCancer-relevant:              {df['cancer_relevant'].sum()}")
    print(f"\nChromosome distribution (top 15):")
    print(df["chrom"].value_counts().head(15).to_string())

    if (df["gene_type"] == "protein_coding").sum() > 0:
        print(f"\nProtein-coding novel genes:")
        coding = df[df["gene_type"] == "protein_coding"].sort_values(["chrom","start"])
        print(coding[["chrom","start","end","name","novel_type",
                       "extra_copy_num","cancer_relevant"]].to_string(index=False))

    # ── Save BED ─────────────────────────────────────────────
    bed = df[["chrom","start","end","name"]].copy()
    bed["score"]  = df["cancer_relevant"].astype(int)
    bed["strand"] = df["strand"]
    bed.to_csv(args.out, sep="\t", header=False, index=False)
    print(f"\nBED saved:  {args.out}  ({len(bed)} genes)")

    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Full TSV:   {args.out_tsv}")
    print(f"\nNext step:")
    print(f"  Use {args.out} as --dark_bed in 06_intersect_annotate.py")
    print(f"  Use {args.out} as -b in samtools depth")


if __name__ == "__main__":
    main()
