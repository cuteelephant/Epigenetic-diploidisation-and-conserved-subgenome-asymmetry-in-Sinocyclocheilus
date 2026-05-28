#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+", required=True,
                   help="Whitespace-separated tissue=path entries.")
    p.add_argument("--out", required=True, help="Output TSV (gene x tissue TPM).")
    p.add_argument("--collapse", choices=["max", "sum", "mean"], default="max",
                   help="Aggregation across duplicate gene IDs (default: max).")
    return p.parse_args()


def load_stringtie(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    # Be tolerant of header variations
    gene_col = next(c for c in df.columns if c.lower().startswith("gene id"))
    tpm_col = next(c for c in df.columns if c.upper() == "TPM")
    s = df[[gene_col, tpm_col]].rename(
        columns={gene_col: "gene_id", tpm_col: "tpm"}
    )
    return s.set_index("gene_id")["tpm"]


def main():
    args = parse_args()
    series = {}
    for entry in args.inputs:
        if "=" not in entry:
            sys.exit(f"--inputs entries must be tissue=path; got: {entry}")
        tissue, path = entry.split("=", 1)
        s = load_stringtie(Path(path))
        if not s.index.is_unique:
            if args.collapse == "max":
                s = s.groupby(s.index).max()
            elif args.collapse == "sum":
                s = s.groupby(s.index).sum()
            else:
                s = s.groupby(s.index).mean()
        series[tissue] = s
    mat = pd.concat(series, axis=1).fillna(0.0)
    mat.index.name = "gene_id"
    mat.to_csv(args.out, sep="\t")
    print(f"[aggregate_tpm_matrix] genes={mat.shape[0]} tissues={mat.shape[1]} "
          f"-> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
