#!/usr/bin/env python3


from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def load_cov(path: Path) -> pd.DataFrame:
    rows = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rows.append((f[0], int(f[1]), int(f[4]), int(f[4]) + int(f[5])))
    return pd.DataFrame(rows, columns=["chrom", "pos", "meth", "total"])


def load_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None,
                       names=["chrom", "start", "end", "gene", "score", "strand"])


def region_interval(row, region: str, flank: int):
    if region == "promoter":
        if row["strand"] == "+":
            return (row["start"] - flank, row["start"])
        return (row["end"], row["end"] + flank)
    elif region == "body":
        return (row["start"], row["end"])
    raise ValueError(f"Unknown region: {region}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cov",      required=True, type=Path)
    p.add_argument("--gene-bed", required=True, type=Path)
    p.add_argument("--region",   required=True, choices=["promoter", "body"])
    p.add_argument("--flank",    type=int, default=2000)
    p.add_argument("--context",  default="CpG")
    p.add_argument("--outdir",   required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    cov  = load_cov(args.cov)
    bed  = load_bed(args.gene_bed)

    cov_by_chr = {c: g.reset_index(drop=True) for c, g in cov.groupby("chrom")}
    rows = []
    for _, row in bed.iterrows():
        chrom = row["chrom"]
        if chrom not in cov_by_chr:
            rows.append((row["gene"], args.region, None, 0))
            continue
        lo, hi = region_interval(row, args.region, args.flank)
        if hi <= lo:
            rows.append((row["gene"], args.region, None, 0))
            continue
        sub = cov_by_chr[chrom]
        mask = (sub["pos"] >= lo) & (sub["pos"] < hi)
        if not mask.any():
            rows.append((row["gene"], args.region, None, 0))
            continue
        total = int(sub.loc[mask, "total"].sum())
        meth  = int(sub.loc[mask, "meth"].sum())
        level = meth / total if total > 0 else None
        rows.append((row["gene"], args.region, level, total))

    out_df = pd.DataFrame(rows,
                          columns=["gene", "region", "meth_level", "n_reads"])
    out_path = args.outdir / f"gene_{args.region}_methylation.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] wrote {out_path}  ({len(out_df)} genes)")


if __name__ == "__main__":
    main()
