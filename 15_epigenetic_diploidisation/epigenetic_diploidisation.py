#!/usr/bin/env python3


from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Thresholds (matching Methods, "Quantifying epigenetic diploidisation")
# ---------------------------------------------------------------------------
DELTA_METH_THR = 0.3        # |ΔMeth| >= 0.3
LOG2FC_THR     = 1.0        # |log2FC|  >= 1.0
MIN_PROMOTER_READS = 10     # >=10 WGBS reads per promoter
MIN_MAX_TPM = 1.0           # at least one homoeolog with TPM>=1
PSEUDOCOUNT = 1.0           # for log2 ratio

CATEGORIES = [
    "S_silenced",
    "D_silenced",
    "Discordant",
    "Meth_asym_only",
    "Expr_asym_only",
    "Balanced",
]


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------
def classify_pair(delta_meth: float, log2fc: float) -> str:
    """Assign a homoeologous pair to one of six categories.

    Definitions follow Methods section "Quantifying epigenetic diploidisation".
    """
    meth_asym = abs(delta_meth) >= DELTA_METH_THR
    expr_asym = abs(log2fc)     >= LOG2FC_THR

    if meth_asym and expr_asym:
        # Sign-discordant means the MORE methylated copy is the LESS expressed
        # one (consistent with methylation-mediated silencing).
        if delta_meth < 0 and log2fc > 0:
            return "S_silenced"           # SubS hypermethylated & lower expr
        if delta_meth > 0 and log2fc < 0:
            return "D_silenced"           # SubD hypermethylated & lower expr
        return "Discordant"               # same-sign: methylation does NOT
                                          # predict the lower-expressed copy
    if meth_asym and not expr_asym:
        return "Meth_asym_only"
    if expr_asym and not meth_asym:
        return "Expr_asym_only"
    return "Balanced"


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def load_methylation(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"gene", "meth_level", "n_reads"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] {path}: missing columns {missing}")
    df = df[["gene", "meth_level", "n_reads"]].dropna(subset=["gene"])
    df["gene"] = df["gene"].astype(str)
    return df


def load_tpm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"gene", "tpm"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] {path}: missing columns {missing}")
    df = df[["gene", "tpm"]].dropna(subset=["gene"])
    df["gene"] = df["gene"].astype(str)
    return df


def load_pairs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"pair_id", "gene_D", "gene_S"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] {path}: missing columns {missing}")
    df["gene_D"] = df["gene_D"].astype(str)
    df["gene_S"] = df["gene_S"].astype(str)
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def build_table(pairs: pd.DataFrame,
                meth_d: pd.DataFrame, meth_s: pd.DataFrame,
                tpm_d: pd.DataFrame,  tpm_s: pd.DataFrame) -> pd.DataFrame:
    """Join methylation and expression data onto syntenic pairs."""
    df = pairs.copy()
    df = df.merge(meth_d.rename(columns={"gene": "gene_D",
                                         "meth_level": "meth_D",
                                         "n_reads": "nreads_D"}),
                  on="gene_D", how="left")
    df = df.merge(meth_s.rename(columns={"gene": "gene_S",
                                         "meth_level": "meth_S",
                                         "n_reads": "nreads_S"}),
                  on="gene_S", how="left")
    df = df.merge(tpm_d.rename(columns={"gene": "gene_D", "tpm": "tpm_D"}),
                  on="gene_D", how="left")
    df = df.merge(tpm_s.rename(columns={"gene": "gene_S", "tpm": "tpm_S"}),
                  on="gene_S", how="left")
    return df


def apply_qc(df: pd.DataFrame) -> pd.DataFrame:
    """Apply per-pair quality filters."""
    qc = (
        df["meth_D"].notna()
        & df["meth_S"].notna()
        & df["tpm_D"].notna()
        & df["tpm_S"].notna()
        & (df["nreads_D"] >= MIN_PROMOTER_READS)
        & (df["nreads_S"] >= MIN_PROMOTER_READS)
        & (df[["tpm_D", "tpm_S"]].max(axis=1) >= MIN_MAX_TPM)
    )
    return df.loc[qc].copy()


def compute_metrics(df: pd.DataFrame) -> pd.DataFrame:
    df["delta_meth"] = df["meth_D"] - df["meth_S"]
    df["log2fc"]     = np.log2((df["tpm_D"] + PSEUDOCOUNT) /
                               (df["tpm_S"] + PSEUDOCOUNT))
    df["category"]   = [classify_pair(d, l)
                        for d, l in zip(df["delta_meth"], df["log2fc"])]
    return df


def summarise(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    counts = df["category"].value_counts().reindex(CATEGORIES, fill_value=0)
    total = int(counts.sum())
    summary = pd.DataFrame({
        "category": counts.index,
        "count":    counts.values,
        "fraction": counts.values / total if total else 0.0,
    })

    # Spearman correlation across all QC-passed pairs
    rho, p_corr = stats.spearmanr(df["delta_meth"], df["log2fc"])

    # Binomial test on D_silenced vs S_silenced
    n_d = int(counts.get("D_silenced", 0))
    n_s = int(counts.get("S_silenced", 0))
    n_concordant = n_d + n_s
    if n_concordant > 0:
        # scipy >=1.7 returns BinomTestResult; older versions return float
        try:
            binom = stats.binomtest(n_s, n_concordant, p=0.5,
                                    alternative="two-sided")
            p_binom = binom.pvalue
        except AttributeError:                       # very old SciPy
            p_binom = stats.binom_test(n_s, n_concordant, p=0.5,
                                       alternative="two-sided")
        prop_s_silenced = n_s / n_concordant
    else:
        p_binom = float("nan")
        prop_s_silenced = float("nan")

    stats_dict = {
        "n_pairs_QC":              int(total),
        "spearman_rho":            float(rho),
        "spearman_p":              float(p_corr),
        "n_D_silenced":            n_d,
        "n_S_silenced":            n_s,
        "n_concordant":            int(n_concordant),
        "prop_S_silenced":         (float(prop_s_silenced)
                                    if not math.isnan(prop_s_silenced)
                                    else None),
        "binomial_p_twosided":     (float(p_binom)
                                    if not math.isnan(p_binom) else None),
        "thresholds": {
            "delta_meth": DELTA_METH_THR,
            "log2fc":     LOG2FC_THR,
            "min_promoter_reads": MIN_PROMOTER_READS,
            "min_max_tpm": MIN_MAX_TPM,
        },
    }
    return summary, stats_dict


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs",   required=True, type=Path)
    p.add_argument("--meth-d",  required=True, type=Path)
    p.add_argument("--meth-s",  required=True, type=Path)
    p.add_argument("--tpm-d",   required=True, type=Path)
    p.add_argument("--tpm-s",   required=True, type=Path)
    p.add_argument("--species", required=True)
    p.add_argument("--tissue",  required=True)
    p.add_argument("--outdir",  required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    pairs  = load_pairs(args.pairs)
    meth_d = load_methylation(args.meth_d)
    meth_s = load_methylation(args.meth_s)
    tpm_d  = load_tpm(args.tpm_d)
    tpm_s  = load_tpm(args.tpm_s)

    raw = build_table(pairs, meth_d, meth_s, tpm_d, tpm_s)
    qc  = apply_qc(raw)
    qc  = compute_metrics(qc)

    summary, stats_dict = summarise(qc)

    stem = f"{args.species}_{args.tissue}"
    qc.to_csv(args.outdir / f"{stem}_pairs_classified.tsv",
              sep="\t", index=False)
    summary.to_csv(args.outdir / f"{stem}_summary.tsv",
                   sep="\t", index=False)
    with open(args.outdir / f"{stem}_stats.json", "w") as fh:
        json.dump(stats_dict, fh, indent=2)

    # Console report mirrors numbers given in the paper
    print(f"[{args.species} | {args.tissue}] "
          f"QC-passed pairs: {stats_dict['n_pairs_QC']:,}")
    print(f"  Spearman rho = {stats_dict['spearman_rho']:.4f}, "
          f"P = {stats_dict['spearman_p']:.2e}")
    print(f"  S_silenced = {stats_dict['n_S_silenced']}, "
          f"D_silenced = {stats_dict['n_D_silenced']}, "
          f"binomial P = {stats_dict['binomial_p_twosided']:.2e}")


if __name__ == "__main__":
    main()
