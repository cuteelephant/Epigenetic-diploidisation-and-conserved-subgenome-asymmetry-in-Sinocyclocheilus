#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SUPERFAMILIES = ["DNA", "LTR", "LINE", "SINE", "RC", "Retroposon", "Unknown"]
EXCLUDED_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite"}


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_genes(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand"])
    return df.set_index("gene")


def load_te(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"chrom", "start", "end", "te_superfamily", "K"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] {path}: missing columns {missing}")
    df = df[~df["te_superfamily"].isin(EXCLUDED_CLASSES)].copy()
    df["midpoint"] = (df["start"] + df["end"]) // 2
    df["te_superfamily"] = df["te_superfamily"].where(
        df["te_superfamily"].isin(SUPERFAMILIES), "Unknown")
    return df


# ---------------------------------------------------------------------------
# Flanking TE retrieval
# ---------------------------------------------------------------------------
def flanking_te_K(gene_row, te_by_chr: dict[str, pd.DataFrame],
                  flank: int) -> pd.DataFrame:
    """Return TEs whose midpoint is within FLANK of the gene boundaries."""
    chrom = gene_row["chrom"]
    if chrom not in te_by_chr:
        return pd.DataFrame(columns=["te_superfamily", "K"])
    g_start, g_end = int(gene_row["start"]), int(gene_row["end"])
    sub = te_by_chr[chrom]
    mask = (
        ((sub["midpoint"] >= g_start - flank) & (sub["midpoint"] <  g_start)) |
        ((sub["midpoint"] >  g_end)           & (sub["midpoint"] <= g_end + flank))
    )
    return sub.loc[mask, ["te_superfamily", "K"]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bias-table", required=True, type=Path)
    p.add_argument("--gene-bed",   required=True, type=Path)
    p.add_argument("--te-table",   required=True, type=Path)
    p.add_argument("--flank",      type=int, default=2000)
    p.add_argument("--mu",         type=float, default=3.51e-9)
    p.add_argument("--outdir",     required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    bias = pd.read_csv(args.bias_table, sep="\t")
    bias = bias[bias["category"].isin(["SubD_biased", "SubS_biased"])].copy()
    genes = load_genes(args.gene_bed)
    te    = load_te(args.te_table)
    te_by_chr = {c: g.reset_index(drop=True) for c, g in te.groupby("chrom")}

    # --- per-pair medians + pooled K records -----------------------------
    per_pair = []
    pooled_rows = []
    for _, row in bias.iterrows():
        gd, gs = row["gene_D"], row["gene_S"]
        if gd not in genes.index or gs not in genes.index:
            continue
        te_d = flanking_te_K(genes.loc[gd], te_by_chr, args.flank)
        te_s = flanking_te_K(genes.loc[gs], te_by_chr, args.flank)
        # assign role
        if row["category"] == "SubD_biased":
            dom_te, sub_te = te_d, te_s
        else:
            dom_te, sub_te = te_s, te_d

        if dom_te.empty and sub_te.empty:
            continue

        per_pair.append({
            "pair_id":        row["pair_id"],
            "category":       row["category"],
            "median_K_dom":   float(np.median(dom_te["K"])) if not dom_te.empty else np.nan,
            "median_K_sub":   float(np.median(sub_te["K"])) if not sub_te.empty else np.nan,
            "n_te_dom":       int(len(dom_te)),
            "n_te_sub":       int(len(sub_te)),
        })

        for _, tk in dom_te.iterrows():
            pooled_rows.append({"role": "dominant",
                                "te_superfamily": tk["te_superfamily"],
                                "K": tk["K"]})
        for _, tk in sub_te.iterrows():
            pooled_rows.append({"role": "submissive",
                                "te_superfamily": tk["te_superfamily"],
                                "K": tk["K"]})

    per_pair_df = pd.DataFrame(per_pair)
    pooled_df   = pd.DataFrame(pooled_rows)
    per_pair_df.to_csv(args.outdir / "te_age_per_pair.tsv",
                       sep="\t", index=False)
    pooled_df.to_csv(args.outdir / "te_age_pooled.tsv",
                     sep="\t", index=False)

    # --- paired Wilcoxon across pairs (all superfamilies pooled) ---------
    paired_input = per_pair_df.dropna(subset=["median_K_dom", "median_K_sub"])
    if not paired_input.empty:
        stat, p_val = stats.wilcoxon(paired_input["median_K_sub"],
                                     paired_input["median_K_dom"],
                                     alternative="less")
    else:
        stat, p_val = float("nan"), float("nan")

    # --- per-superfamily Mann-Whitney (unpaired) on pooled K -------------
    sf_stats = {}
    for sf in SUPERFAMILIES:
        d_k = pooled_df.loc[(pooled_df["role"] == "dominant") &
                            (pooled_df["te_superfamily"] == sf), "K"].values
        s_k = pooled_df.loc[(pooled_df["role"] == "submissive") &
                            (pooled_df["te_superfamily"] == sf), "K"].values
        if len(d_k) < 5 or len(s_k) < 5:
            continue
        u, p = stats.mannwhitneyu(s_k, d_k, alternative="less")
        delta_k = float(np.median(d_k) - np.median(s_k))
        delta_t_myr = delta_k / (2 * args.mu) / 1e6
        sf_stats[sf] = {
            "n_dominant":     int(len(d_k)),
            "n_submissive":   int(len(s_k)),
            "median_K_dom":   float(np.median(d_k)),
            "median_K_sub":   float(np.median(s_k)),
            "delta_K":        delta_k,
            "delta_T_Myr":    delta_t_myr,
            "U":              float(u),
            "p_less":         float(p),
        }

    out_stats = {
        "paired_wilcoxon": {
            "n_pairs":  int(len(paired_input)),
            "W":        float(stat) if not np.isnan(stat) else None,
            "p_less":   float(p_val) if not np.isnan(p_val) else None,
        },
        "mu_subs_site_year": args.mu,
        "per_superfamily":    sf_stats,
    }
    with open(args.outdir / "te_age_stats.json", "w") as fh:
        json.dump(out_stats, fh, indent=2)
    print(json.dumps(out_stats, indent=2))


if __name__ == "__main__":
    main()
