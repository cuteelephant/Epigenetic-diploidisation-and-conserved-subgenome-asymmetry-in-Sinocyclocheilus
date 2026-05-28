#!/usr/bin/env python3


from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO          # preferred
    HAVE_BIOPYTHON = True
except Exception:
    HAVE_BIOPYTHON = False


def parse_fasta(path: Path):
    """Minimal fasta parser yielding (header, seq). Avoids hard biopython dep."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq)


def build_cds_index(cds_dir: Path) -> dict[str, str]:
    index: dict[str, str] = {}
    for fa in cds_dir.glob("*"):
        if fa.suffix.lower() not in {".fa", ".fasta", ".fna", ".cds"}:
            continue
        for name, seq in parse_fasta(fa):
            index[name] = seq
    return index


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--aln-fa",  required=True, type=Path)
    p.add_argument("--cds-dir", required=True, type=Path)
    p.add_argument("--out",     required=True, type=Path)
    args = p.parse_args()

    index = build_cds_index(args.cds_dir)
    written = 0
    with open(args.out, "w") as oh:
        for name, _ in parse_fasta(args.aln_fa):
            cds = index.get(name)
            if cds is None:
                print(f"[WARN] missing CDS for {name}", file=sys.stderr)
                continue
            oh.write(f">{name}\n{cds}\n")
            written += 1
    print(f"[INFO] wrote {written} CDS records to {args.out}")


if __name__ == "__main__":
    main()
