#!/usr/bin/env python3


from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# Plot is optional — fall back to "no plot" if matplotlib missing
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_coverage(path: Path) -> pd.DataFrame:
    """Bismark .cov format: chrom start end pct_meth count_meth count_unmeth."""
    rows = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            chrom = f[0]
            pos   = int(f[1])         # 1-based in Bismark cov
            cm    = int(f[4])
            cu    = int(f[5])
            rows.append((chrom, pos, cm, cm + cu))
    df = pd.DataFrame(rows, columns=["chrom", "pos", "meth", "total"])
    return df


def load_features(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "name", "score", "strand"])
    return df


def load_subgenome_map(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["name", "subgenome"])
    return dict(zip(df["name"], df["subgenome"]))


# ---------------------------------------------------------------------------
# Weighted methylation per bin per feature
# ---------------------------------------------------------------------------
def feature_profile(feature, cov_chr_pos, flank: int, n_bins: int) -> np.ndarray:
    """Return length 3*n_bins vector of weighted methylation."""
    chrom  = feature["chrom"]
    start  = int(feature["start"])
    end    = int(feature["end"])
    strand = feature["strand"]
    if chrom not in cov_chr_pos:
        return np.full(3 * n_bins, np.nan)

    df = cov_chr_pos[chrom]
    pos    = df["pos"].values
    meth   = df["meth"].values
    total  = df["total"].values

    profile = np.full(3 * n_bins, np.nan)

    # define three intervals: upstream, body, downstream (TSS/TTS strand-aware)
    body_len = end - start
    if body_len <= 0:
        return profile
    bins_def = [
        (start - flank, start,    "up"),
        (start,         end,      "body"),
        (end,           end+flank, "down"),
    ]

    for k, (lo, hi, _) in enumerate(bins_def):
        if hi <= lo:
            continue
        idx = (pos >= lo) & (pos < hi)
        if not idx.any():
            continue
        p_sel = pos[idx]
        m_sel = meth[idx]
        t_sel = total[idx]
        # uniform binning
        edges = np.linspace(lo, hi, n_bins + 1)
        # assign each cytosine to a bin
        bin_id = np.minimum(np.searchsorted(edges, p_sel, side="right") - 1,
                            n_bins - 1)
        bin_meth  = np.bincount(bin_id, weights=m_sel, minlength=n_bins)
        bin_total = np.bincount(bin_id, weights=t_sel, minlength=n_bins)
        vals = np.full(n_bins, np.nan)
        nz   = bin_total > 0
        vals[nz] = bin_meth[nz] / bin_total[nz]
        # write into profile slot
        profile[k * n_bins:(k + 1) * n_bins] = vals

    # reverse for minus strand so all features are 5' -> 3'
    if strand == "-":
        # reverse each region within its own slot
        for k in range(3):
            seg = profile[k * n_bins:(k + 1) * n_bins][::-1]
            profile[k * n_bins:(k + 1) * n_bins] = seg
        # also flip up/down order: what is "downstream" on minus is "upstream"
        # of the gene. swap region 0 and 2.
        up   = profile[0:n_bins].copy()
        down = profile[2 * n_bins:3 * n_bins].copy()
        profile[0:n_bins]            = down
        profile[2 * n_bins:3 * n_bins] = up
    return profile


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cov",           required=True, type=Path)
    p.add_argument("--features-bed",  required=True, type=Path)
    p.add_argument("--subgenome-map", required=True, type=Path)
    p.add_argument("--flank",         type=int, default=2000)
    p.add_argument("--bins",          type=int, default=20)
    p.add_argument("--outdir",        required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    cov = load_coverage(args.cov)
    cov_by_chr = {c: g.reset_index(drop=True) for c, g in cov.groupby("chrom")}
    feats = load_features(args.features_bed)
    sub_map = load_subgenome_map(args.subgenome_map)

    n_bins = args.bins
    matrix = []
    name_keep = []
    sub_keep = []
    for _, feat in feats.iterrows():
        if feat["name"] not in sub_map:
            continue
        prof = feature_profile(feat, cov_by_chr, args.flank, n_bins)
        if np.all(np.isnan(prof)):
            continue
        matrix.append(prof)
        name_keep.append(feat["name"])
        sub_keep.append(sub_map[feat["name"]])

    if not matrix:
        sys.exit("[ERROR] no features had coverage data")

    M = np.asarray(matrix)
    cols = ([f"up_{i}"   for i in range(n_bins)] +
            [f"body_{i}" for i in range(n_bins)] +
            [f"down_{i}" for i in range(n_bins)])
    df_mat = pd.DataFrame(M, columns=cols)
    df_mat.insert(0, "subgenome", sub_keep)
    df_mat.insert(0, "name",      name_keep)
    df_mat.to_csv(args.outdir / "metaplot_matrix.tsv",
                  sep="\t", index=False)

    # per-bin mean per subgenome
    mean_d = df_mat.loc[df_mat["subgenome"] == "SubD", cols].mean(axis=0)
    mean_s = df_mat.loc[df_mat["subgenome"] == "SubS", cols].mean(axis=0)
    pd.DataFrame({"bin": cols,
                  "mean_SubD": mean_d.values,
                  "mean_SubS": mean_s.values}).to_csv(
        args.outdir / "metaplot_mean.tsv", sep="\t", index=False)

    # per-bin Wilcoxon rank-sum (Mann-Whitney U) SubD vs SubS
    pvals = []
    for c in cols:
        d_vals = df_mat.loc[df_mat["subgenome"] == "SubD", c].dropna()
        s_vals = df_mat.loc[df_mat["subgenome"] == "SubS", c].dropna()
        if len(d_vals) < 5 or len(s_vals) < 5:
            pvals.append(np.nan)
            continue
        _, p_val = stats.mannwhitneyu(d_vals, s_vals, alternative="two-sided")
        pvals.append(p_val)
    pd.DataFrame({"bin": cols, "p_value": pvals}).to_csv(
        args.outdir / "metaplot_pvals.tsv", sep="\t", index=False)

    if HAS_MPL:
        x = np.arange(3 * n_bins)
        plt.figure(figsize=(7, 3.5))
        plt.plot(x, mean_d.values, label="SubD", linewidth=2)
        plt.plot(x, mean_s.values, label="SubS", linewidth=2)
        plt.axvline(n_bins,     color="grey", linestyle=":")
        plt.axvline(2 * n_bins, color="grey", linestyle=":")
        plt.xticks([0, n_bins, 2 * n_bins, 3 * n_bins],
                   ["-2kb", "TSS", "TTS", "+2kb"])
        plt.ylabel("Weighted CG methylation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.outdir / "metaplot.pdf")
    else:
        print("[INFO] matplotlib not available; skipping PDF.")


if __name__ == "__main__":
    main()
