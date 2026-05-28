#!/usr/bin/env python3


from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

try:
    from sklearn.mixture import GaussianMixture
    HAVE_GMM = True
except Exception:
    HAVE_GMM = False

try:
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


def load_ks(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    ks_col = next((c for c in df.columns if c.lower().startswith("ks")), None)
    if ks_col is None:
        sys.exit("[ERROR] no Ks column found")
    df = df.rename(columns={ks_col: "Ks"})
    df["Ks"] = pd.to_numeric(df["Ks"], errors="coerce")
    return df.dropna(subset=["Ks"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ks",           required=True, type=Path, nargs="+")
    p.add_argument("--labels",       nargs="+", default=None,
                   help="one label per --ks input")
    p.add_argument("--out",          required=True, type=Path)
    p.add_argument("--max-ks",       type=float, default=5.0)
    p.add_argument("--n-components", type=int,   default=4)
    args = p.parse_args()

    if args.labels and len(args.labels) != len(args.ks):
        sys.exit("[ERROR] --labels must match --ks count")

    frames = []
    for i, ks_path in enumerate(args.ks):
        df = load_ks(ks_path)
        df["comparison"] = args.labels[i] if args.labels else ks_path.stem
        frames.append(df[["Ks", "comparison"]])
    ks_df = pd.concat(frames, ignore_index=True)
    ks_df = ks_df[(ks_df["Ks"] > 0) & (ks_df["Ks"] < args.max_ks)]

    # KDE per comparison
    grid = np.linspace(0, args.max_ks, 1000)
    density_rows = []
    for label, g in ks_df.groupby("comparison"):
        if len(g) < 30:
            continue
        kde = gaussian_kde(g["Ks"].values, bw_method=0.15)
        dens = kde(grid)
        for x, y in zip(grid, dens):
            density_rows.append({"comparison": label, "Ks": x, "density": y})
    density_df = pd.DataFrame(density_rows)
    density_df.to_csv(f"{args.out}.density.tsv", sep="\t", index=False)

    # Gaussian-mixture peak identification on pooled distribution
    peak_rows = []
    if HAVE_GMM:
        x = ks_df["Ks"].values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=args.n_components,
                              covariance_type="full", random_state=0)
        gmm.fit(x)
        means = gmm.means_.flatten()
        sds   = np.sqrt(gmm.covariances_.flatten())
        wts   = gmm.weights_
        order = np.argsort(means)
        for k in order:
            peak_rows.append({"component": int(k),
                              "mean_Ks":   float(means[k]),
                              "sd_Ks":     float(sds[k]),
                              "weight":    float(wts[k])})
    pd.DataFrame(peak_rows).to_csv(f"{args.out}.peaks.tsv",
                                   sep="\t", index=False)

    if HAVE_MPL:
        plt.figure(figsize=(7, 4))
        for label, g in density_df.groupby("comparison"):
            plt.plot(g["Ks"], g["density"], label=label)
        # Annotate the manuscript-reported peaks
        for peak_x, name in [(0.05, "recent"),
                             (0.13, "CS4R"),
                             (2.00, "TS3R")]:
            if peak_x <= args.max_ks:
                plt.axvline(peak_x, color="grey", linestyle=":", linewidth=0.8)
                plt.text(peak_x, plt.ylim()[1] * 0.9, name,
                         ha="center", fontsize=8, color="grey")
        plt.xlabel("Ks")
        plt.ylabel("Density")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(f"{args.out}.pdf")

    print(f"[INFO] density -> {args.out}.density.tsv")
    print(f"[INFO] peaks   -> {args.out}.peaks.tsv")


if __name__ == "__main__":
    main()
