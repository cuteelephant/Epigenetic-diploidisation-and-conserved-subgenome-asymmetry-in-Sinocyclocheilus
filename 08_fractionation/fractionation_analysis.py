#!/usr/bin/env python3


from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quota-output", required=True, type=Path)
    p.add_argument("--ref-bed",      required=True, type=Path)
    p.add_argument("--window-genes", type=int, default=100)
    p.add_argument("--out",          required=True, type=Path)
    args = p.parse_args()

    qa = pd.read_csv(args.quota_output, sep="\t",
                     names=["ref_gene", "query_gene", "subgenome"])
    ref = pd.read_csv(args.ref_bed, sep="\t", header=None,
                      names=["chrom", "start", "end", "gene", "score", "strand"])
    # ordered by chromosome and start
    ref = ref.sort_values(["chrom", "start"]).reset_index(drop=True)
    ref["order"] = range(len(ref))

    # presence flags
    q_d = qa[qa["subgenome"] == "SubD"]["ref_gene"].unique()
    q_s = qa[qa["subgenome"] == "SubS"]["ref_gene"].unique()
    ref["has_D"] = ref["gene"].isin(q_d).astype(int)
    ref["has_S"] = ref["gene"].isin(q_s).astype(int)

    # sliding window per chromosome
    rows = []
    w = args.window_genes
    for chrom, g in ref.groupby("chrom"):
        g = g.reset_index(drop=True)
        for i in range(0, len(g) - w + 1, max(1, w // 2)):
            window = g.iloc[i:i + w]
            rows.append({
                "chrom": chrom,
                "window_start_gene": window["gene"].iloc[0],
                "window_end_gene":   window["gene"].iloc[-1],
                "ret_SubD":  window["has_D"].mean(),
                "ret_SubS":  window["has_S"].mean(),
            })
    win_df = pd.DataFrame(rows)
    win_df.to_csv(f"{args.out}.windows.tsv", sep="\t", index=False)

    overall = {
        "n_ref_genes":        int(len(ref)),
        "retention_SubD":     float(ref["has_D"].mean()),
        "retention_SubS":     float(ref["has_S"].mean()),
        "delta_D_minus_S":    float(ref["has_D"].mean() - ref["has_S"].mean()),
    }
    pd.DataFrame([overall]).to_csv(f"{args.out}.summary.tsv",
                                   sep="\t", index=False)
    print(pd.DataFrame([overall]).to_string(index=False))


if __name__ == "__main__":
    main()
