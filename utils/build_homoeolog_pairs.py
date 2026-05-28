#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def load_subgenome(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["gene_id", "subgenome"], comment="#")
    return dict(zip(df.gene_id, df.subgenome.str.upper()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--subgenome", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sg = load_subgenome(args.subgenome)
    rows = []
    block = 0
    dropped = 0
    with open(args.anchors) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("###"):
                block += 1
                continue
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            g1, g2 = parts[0], parts[1]
            score = parts[2] if len(parts) > 2 else "NA"
            sg1, sg2 = sg.get(g1), sg.get(g2)
            if sg1 == "D" and sg2 == "S":
                gene_D, gene_S = g1, g2
            elif sg1 == "S" and sg2 == "D":
                gene_D, gene_S = g2, g1
            else:
                dropped += 1
                continue
            pair_id = f"{gene_D}__{gene_S}"
            rows.append((pair_id, gene_D, gene_S, f"block{block:04d}", score))

    df = pd.DataFrame(rows, columns=["pair_id", "gene_D", "gene_S",
                                     "block_id", "score"])
    # Keep one row per (gene_D, gene_S) — anchors files can repeat
    df = df.drop_duplicates(subset=["gene_D", "gene_S"])
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[build_homoeolog_pairs] pairs kept = {len(df)} "
          f"(dropped {dropped} mis-labelled / cross-subgenome)", file=sys.stderr)


if __name__ == "__main__":
    main()
