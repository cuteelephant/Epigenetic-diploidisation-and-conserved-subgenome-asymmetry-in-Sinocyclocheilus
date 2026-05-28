#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True,
                    help="TSV with one row per feature, methylation columns.")
    ap.add_argument("--group-col", default=None,
                    help="Column that classifies rows (e.g. subgenome).")
    ap.add_argument("--value-cols", nargs="+",
                    help="Columns to use as methylation values "
                         "(default: all numeric columns).")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--cmap", default="viridis")
    args = ap.parse_args()

    df = pd.read_csv(args.matrix, sep="\t")
    if args.value_cols:
        value_cols = args.value_cols
    else:
        value_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if args.group_col and args.group_col in value_cols:
            value_cols.remove(args.group_col)
    if not value_cols:
        sys.exit("No numeric methylation columns found.")

    mat = df[value_cols].copy()

    # Heatmap (subsample to 5000 rows for tractable clustering)
    sub = mat.sample(min(5000, len(mat)), random_state=42)
    g = sns.clustermap(sub.fillna(0), cmap=args.cmap, vmin=0, vmax=1,
                       figsize=(max(4, len(value_cols) * 0.6), 8),
                       yticklabels=False)
    g.fig.suptitle("Promoter methylation (subsampled)", y=1.02)
    g.savefig(f"{args.out_prefix}_heatmap.pdf", bbox_inches="tight")
    plt.close("all")

    # Group-stratified ridge / violin
    if args.group_col and args.group_col in df.columns:
        long = df.melt(id_vars=[args.group_col], value_vars=value_cols,
                       var_name="sample", value_name="meth")
        fig, ax = plt.subplots(figsize=(max(4, len(value_cols) * 0.7), 4))
        sns.violinplot(data=long, x="sample", y="meth",
                       hue=args.group_col, split=True, inner="quartile",
                       ax=ax, cut=0)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Methylation level")
        ax.set_xlabel("")
        plt.xticks(rotation=30, ha="right")
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}_violin.pdf")
        plt.close(fig)

        # Stats summary
        rows = []
        for col in value_cols:
            groups = [g[col].dropna().values for _, g in df.groupby(args.group_col)]
            if len(groups) == 2 and all(len(x) > 0 for x in groups):
                u, p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
                rows.append({
                    "sample": col,
                    "n1": len(groups[0]), "n2": len(groups[1]),
                    "mean1": np.mean(groups[0]), "mean2": np.mean(groups[1]),
                    "U": u, "p": p,
                })
        if rows:
            pd.DataFrame(rows).to_csv(f"{args.out_prefix}_stats.tsv",
                                      sep="\t", index=False)

    print(f"[plot_methylation_heatmap] wrote {args.out_prefix}_*.pdf",
          file=sys.stderr)


if __name__ == "__main__":
    main()
