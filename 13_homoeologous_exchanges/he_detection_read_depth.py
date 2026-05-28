#!/usr/bin/env python3


from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_depth(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "mean_depth"])
    return df


def load_synteny(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom_a", "start_a", "end_a",
                            "chrom_b", "start_b", "end_b"])
    return df


def flag_outliers(depth: pd.DataFrame,
                  high_range: tuple[float, float],
                  low_range:  tuple[float, float]) -> pd.DataFrame:
    out = []
    for chrom, g in depth.groupby("chrom"):
        # chromosome-wide mean ignoring zero-depth windows (likely gaps)
        nonzero = g.loc[g["mean_depth"] > 0, "mean_depth"]
        if nonzero.empty:
            continue
        chrom_mean = float(np.mean(nonzero))
        g = g.copy()
        g["ratio"] = g["mean_depth"] / chrom_mean
        g["flag"]  = "normal"
        g.loc[g["ratio"].between(*high_range, inclusive="both"), "flag"] = "high"
        g.loc[g["ratio"].between(*low_range,  inclusive="both"), "flag"] = "low"
        g["chrom_mean"] = chrom_mean
        out.append(g)
    return pd.concat(out, ignore_index=True)


def merge_runs(flagged: pd.DataFrame) -> pd.DataFrame:
    """Merge contiguous windows with the same non-normal flag into blocks."""
    blocks = []
    for chrom, g in flagged.groupby("chrom"):
        g = g.sort_values("start").reset_index(drop=True)
        cur_flag = None
        cur_start = cur_end = None
        for _, row in g.iterrows():
            if row["flag"] == "normal":
                if cur_flag is not None:
                    blocks.append((chrom, cur_start, cur_end, cur_flag))
                cur_flag = None
                continue
            if cur_flag is None:
                cur_flag = row["flag"]
                cur_start, cur_end = row["start"], row["end"]
            elif row["flag"] == cur_flag and row["start"] == cur_end:
                cur_end = row["end"]
            else:
                blocks.append((chrom, cur_start, cur_end, cur_flag))
                cur_flag = row["flag"]
                cur_start, cur_end = row["start"], row["end"]
        if cur_flag is not None:
            blocks.append((chrom, cur_start, cur_end, cur_flag))
    return pd.DataFrame(blocks, columns=["chrom", "start", "end", "flag"])


def within_synteny(block: pd.Series, synteny: pd.DataFrame,
                   min_overlap: float = 0.6) -> bool:
    """Block overlaps any syntenic block on side A by >= min_overlap fraction."""
    b_len = block["end"] - block["start"]
    if b_len <= 0:
        return False
    matches = synteny[(synteny["chrom_a"] == block["chrom"]) &
                      (synteny["end_a"]   >  block["start"]) &
                      (synteny["start_a"] <  block["end"])]
    if matches.empty:
        return False
    for _, m in matches.iterrows():
        ov = max(0, min(block["end"], m["end_a"]) - max(block["start"], m["start_a"]))
        if ov / b_len >= min_overlap:
            return True
    return False


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--depth-bed",   required=True, type=Path)
    p.add_argument("--synteny-bed", required=True, type=Path)
    p.add_argument("--min-len",     type=int,   default=60000)
    p.add_argument("--high-range",  default="1.5,5.0")
    p.add_argument("--low-range",   default="0.0,0.5")
    p.add_argument("--outdir",      required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    depth = load_depth(args.depth_bed)
    synteny = load_synteny(args.synteny_bed)

    high_rng = tuple(float(x) for x in args.high_range.split(","))
    low_rng  = tuple(float(x) for x in args.low_range.split(","))

    flagged = flag_outliers(depth, high_rng, low_rng)
    flagged[flagged["flag"] != "normal"].to_csv(
        args.outdir / "HE_candidate_windows.bed",
        sep="\t", index=False)

    blocks = merge_runs(flagged)
    blocks["length"] = blocks["end"] - blocks["start"]
    blocks = blocks[blocks["length"] >= args.min_len].reset_index(drop=True)
    blocks["in_synteny"] = blocks.apply(
        lambda r: within_synteny(r, synteny), axis=1)
    blocks_final = blocks[blocks["in_synteny"]].copy()

    blocks_final[["chrom", "start", "end", "flag", "length"]].to_csv(
        args.outdir / "HE_blocks_filtered.bed",
        sep="\t", index=False, header=False)

    summary = {
        "n_candidate_windows": int((flagged["flag"] != "normal").sum()),
        "n_blocks_before_synteny_filter":
            int(len(blocks)),
        "n_HE_blocks_final":   int(len(blocks_final)),
        "high_range_ratio":    list(high_rng),
        "low_range_ratio":     list(low_rng),
        "min_block_length_bp": args.min_len,
        "per_chrom_counts":    blocks_final["chrom"].value_counts().to_dict(),
    }
    with open(args.outdir / "HE_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
