#!/usr/bin/env python3


from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

FOLD_CHANGE_THR = 2.0           # >2-fold TPM ratio defines bias
MIN_TPM_ANY_SAMPLE = 1.0        # exclude pairs with TPM<1 in all samples
PSEUDOCOUNT = 1.0


def load_pairs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col in ("pair_id", "gene_D", "gene_S"):
        if col not in df.columns:
            sys.exit(f"[ERROR] {path}: missing column '{col}'")
    df["gene_D"] = df["gene_D"].astype(str)
    df["gene_S"] = df["gene_S"].astype(str)
    return df


def load_tpm_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.columns[0] != "gene":
        df = df.rename(columns={df.columns[0]: "gene"})
    df["gene"] = df["gene"].astype(str)
    df = df.set_index("gene")
    return df


def classify_bias(tpm_d: float, tpm_s: float) -> str:
    if tpm_d < MIN_TPM_ANY_SAMPLE and tpm_s < MIN_TPM_ANY_SAMPLE:
        return "low_expr"
    # use pseudocount-protected ratios in both directions
    if (tpm_d + PSEUDOCOUNT) / (tpm_s + PSEUDOCOUNT) > FOLD_CHANGE_THR:
        return "SubD_biased"
    if (tpm_s + PSEUDOCOUNT) / (tpm_d + PSEUDOCOUNT) > FOLD_CHANGE_THR:
        return "SubS_biased"
    return "balanced"


def per_tissue_table(pairs: pd.DataFrame,
                     tpm: pd.DataFrame,
                     tissue: str) -> pd.DataFrame:
    if tissue not in tpm.columns:
        sys.exit(f"[ERROR] tissue '{tissue}' not in TPM matrix")
    rows = []
    for _, row in pairs.iterrows():
        gd, gs = row["gene_D"], row["gene_S"]
        td = float(tpm.at[gd, tissue]) if gd in tpm.index else np.nan
        ts = float(tpm.at[gs, tissue]) if gs in tpm.index else np.nan
        if np.isnan(td) or np.isnan(ts):
            continue
        rows.append({
            "pair_id":  row["pair_id"],
            "gene_D":   gd,
            "gene_S":   gs,
            "tpm_D":    td,
            "tpm_S":    ts,
            "log2fc":   np.log2((td + PSEUDOCOUNT) / (ts + PSEUDOCOUNT)),
            "category": classify_bias(td, ts),
        })
    return pd.DataFrame(rows)


def summarise(per_tissue: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Wilcoxon signed-rank on per-tissue log2FC distribution."""
    out = []
    for tissue, df in per_tissue.items():
        sub = df[df["category"] != "low_expr"]
        if sub.empty:
            continue
        # one-sample Wilcoxon: is log2FC distribution shifted from 0?
        try:
            stat, p = stats.wilcoxon(sub["log2fc"].values,
                                     zero_method="wilcox",
                                     alternative="two-sided")
        except ValueError:
            stat, p = float("nan"), float("nan")
        n_d = int((sub["category"] == "SubD_biased").sum())
        n_s = int((sub["category"] == "SubS_biased").sum())
        n_bal = int((sub["category"] == "balanced").sum())
        out.append({
            "tissue":         tissue,
            "n_informative":  len(sub),
            "n_SubD_biased":  n_d,
            "n_SubS_biased":  n_s,
            "n_balanced":     n_bal,
            "median_log2fc":  float(np.median(sub["log2fc"])),
            "wilcoxon_stat":  float(stat) if not np.isnan(stat) else None,
            "wilcoxon_p":     float(p)    if not np.isnan(p)    else None,
        })
    return pd.DataFrame(out)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs",      required=True, type=Path)
    p.add_argument("--tpm-matrix", required=True, type=Path)
    p.add_argument("--tissues",    required=True, nargs="+",
                   help="tissue column names to analyse")
    p.add_argument("--ks-pairs",   type=Path, default=None,
                   help="optional Ks table; pairs with Ks < --ks-max kept")
    p.add_argument("--ks-max",     type=float, default=0.05)
    p.add_argument("--species",    required=True)
    p.add_argument("--outdir",     required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    pairs = load_pairs(args.pairs)
    if args.ks_pairs is not None:
        ks = pd.read_csv(args.ks_pairs, sep="\t")
        if "pair_id" not in ks.columns or "Ks" not in ks.columns:
            sys.exit("[ERROR] --ks-pairs must have columns pair_id, Ks")
        keep_ids = set(ks.loc[ks["Ks"] < args.ks_max, "pair_id"])
        pairs = pairs[pairs["pair_id"].isin(keep_ids)].reset_index(drop=True)
        print(f"[INFO] Restricted to {len(pairs)} pairs with Ks < {args.ks_max}")

    tpm = load_tpm_matrix(args.tpm_matrix)

    per_tissue: dict[str, pd.DataFrame] = {}
    long_rows = []
    for tissue in args.tissues:
        df = per_tissue_table(pairs, tpm, tissue)
        per_tissue[tissue] = df
        df.to_csv(args.outdir / f"{args.species}_HEB_{tissue}.tsv",
                  sep="\t", index=False)
        df_long = df[["pair_id", "log2fc", "category"]].copy()
        df_long["tissue"] = tissue
        long_rows.append(df_long)

    pd.concat(long_rows, ignore_index=True).to_csv(
        args.outdir / f"{args.species}_HEB_log2fc_matrix.tsv",
        sep="\t", index=False)

    summary = summarise(per_tissue)
    summary.to_csv(args.outdir / f"{args.species}_HEB_summary.tsv",
                   sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
