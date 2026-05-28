#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Coordinate loading
# ---------------------------------------------------------------------------
def load_gene_coords(path: Path) -> pd.DataFrame:
    """Return gene -> (chrom, start, end, strand) DataFrame.

    Supports GFF3 or 6-column BED. Detection is by extension.
    """
    suffix = path.suffix.lower()
    if suffix in {".gff", ".gff3"}:
        rows = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                attrs = dict(item.split("=", 1) for item in parts[8].split(";")
                             if "=" in item)
                gid = attrs.get("ID") or attrs.get("Name")
                if not gid:
                    continue
                rows.append((gid, parts[0], int(parts[3]) - 1,
                             int(parts[4]), parts[6]))
        return pd.DataFrame(rows, columns=["gene", "chrom", "start",
                                            "end", "strand"]).set_index("gene")
    # assume BED-like
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand"])
    return df.set_index("gene")


def load_te(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None,
                     usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "name"])
    return df


# ---------------------------------------------------------------------------
# Flank TE coverage
# ---------------------------------------------------------------------------
def flanking_te_fraction(gene_row, te_by_chr: dict[str, pd.DataFrame],
                         flank: int) -> float:
    chrom = gene_row["chrom"]
    if chrom not in te_by_chr:
        return 0.0
    g_start, g_end = int(gene_row["start"]), int(gene_row["end"])
    up_lo, up_hi = max(0, g_start - flank), g_start
    dn_lo, dn_hi = g_end, g_end + flank

    sub = te_by_chr[chrom]
    # naive overlap: rely on pandas; TE list per chromosome is short enough
    overlap_up = sub[(sub["end"] > up_lo) & (sub["start"] < up_hi)]
    overlap_dn = sub[(sub["end"] > dn_lo) & (sub["start"] < dn_hi)]

    bp_covered = 0
    for _, te in overlap_up.iterrows():
        bp_covered += max(0, min(te["end"], up_hi) - max(te["start"], up_lo))
    for _, te in overlap_dn.iterrows():
        bp_covered += max(0, min(te["end"], dn_hi) - max(te["start"], dn_lo))
    return bp_covered / (2 * flank)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bias-table", required=True, type=Path)
    p.add_argument("--gff-genes",  required=True, type=Path)
    p.add_argument("--bed-te",     required=True, type=Path)
    p.add_argument("--flank",      type=int, default=2000)
    p.add_argument("--outdir",     required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    bias = pd.read_csv(args.bias_table, sep="\t")
    if "category" not in bias.columns:
        sys.exit("[ERROR] --bias-table must contain 'category' column")
    bias = bias[bias["category"].isin(["SubD_biased", "SubS_biased"])].copy()
    if bias.empty:
        sys.exit("[ERROR] no D_biased/S_biased pairs in bias table")

    genes = load_gene_coords(args.gff_genes)
    te    = load_te(args.bed_te)
    te_by_chr = {c: g.reset_index(drop=True) for c, g in te.groupby("chrom")}

    out_rows = []
    for _, row in bias.iterrows():
        gd, gs = row["gene_D"], row["gene_S"]
        if gd not in genes.index or gs not in genes.index:
            continue
        f_d = flanking_te_fraction(genes.loc[gd], te_by_chr, args.flank)
        f_s = flanking_te_fraction(genes.loc[gs], te_by_chr, args.flank)
        # assign dominant / submissive based on bias category
        if row["category"] == "SubD_biased":
            dom, sub = f_d, f_s
        else:
            dom, sub = f_s, f_d
        out_rows.append({
            "pair_id":           row["pair_id"],
            "category":          row["category"],
            "te_frac_dominant":  dom,
            "te_frac_submissive": sub,
        })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.outdir / "flanking_TE_density.tsv",
                  sep="\t", index=False)

    # paired Wilcoxon test (Methods)
    try:
        stat, p_val = stats.wilcoxon(out_df["te_frac_submissive"],
                                     out_df["te_frac_dominant"],
                                     alternative="greater")
    except ValueError:
        stat, p_val = float("nan"), float("nan")

    stats_dict = {
        "n_pairs":          int(len(out_df)),
        "median_dominant":  float(np.median(out_df["te_frac_dominant"])),
        "median_submissive":float(np.median(out_df["te_frac_submissive"])),
        "wilcoxon_W":       float(stat) if not np.isnan(stat) else None,
        "wilcoxon_p_greater":float(p_val) if not np.isnan(p_val) else None,
        "flank_bp":         args.flank,
    }
    with open(args.outdir / "flanking_TE_density_stats.json", "w") as fh:
        json.dump(stats_dict, fh, indent=2)
    print(json.dumps(stats_dict, indent=2))


if __name__ == "__main__":
    main()
