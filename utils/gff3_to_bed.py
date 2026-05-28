#!/usr/bin/env python3

import argparse
import sys
import re
from pathlib import Path


def parse_attrs(s: str) -> dict:
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def load_fai(path: Path) -> dict:
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--feature", default="gene",
                    help="GFF3 feature type, or 'promoter' (default: gene)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--upstream", type=int, default=2000,
                    help="Promoter upstream length when --feature promoter")
    ap.add_argument("--fai", help="genome.fa.fai (needed for promoter clipping)")
    ap.add_argument("--id-key", default="ID",
                    help="Attribute key to use as BED name (default: ID).")
    args = ap.parse_args()

    promoter_mode = (args.feature == "promoter")
    target = "gene" if promoter_mode else args.feature

    sizes = load_fai(Path(args.fai)) if (promoter_mode and args.fai) else {}

    n_in = n_out = 0
    with open(args.gff) as fh, open(args.out, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != target:
                continue
            n_in += 1
            chrom = parts[0]
            start = int(parts[3]) - 1  # GFF3 1-based -> BED 0-based
            end = int(parts[4])
            strand = parts[6] if parts[6] in ("+", "-") else "+"
            attrs = parse_attrs(parts[8])
            name = attrs.get(args.id_key, attrs.get("Name", f"{chrom}:{start}-{end}"))

            if promoter_mode:
                if strand == "+":
                    p_end = start
                    p_start = max(0, p_end - args.upstream)
                else:
                    p_start = end
                    p_end = p_start + args.upstream
                    if chrom in sizes:
                        p_end = min(p_end, sizes[chrom])
                if p_end <= p_start:
                    continue
                out.write(f"{chrom}\t{p_start}\t{p_end}\t{name}\t.\t{strand}\n")
            else:
                out.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")
            n_out += 1

    print(f"[gff3_to_bed] {target}: {n_in} in, {n_out} written -> {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
