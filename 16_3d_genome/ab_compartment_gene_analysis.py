#!/usr/bin/env python3


from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def load_bed(path: Path, names: list[str]) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=names)


def assign_compartment(gene_row, comp_by_chr: dict[str, pd.DataFrame]) -> str:
    chrom = gene_row["chrom"]
    if chrom not in comp_by_chr:
        return "NA"
    centre = (gene_row["start"] + gene_row["end"]) // 2
    sub = comp_by_chr[chrom]
    hits = sub[(sub["start"] <= centre) & (sub["end"] > centre)]
    if hits.empty:
        return "NA"
    return str(hits.iloc[0]["label"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compartments-bed", required=True, type=Path)
    p.add_argument("--gene-bed",         required=True, type=Path)
    p.add_argument("--subgenome-map",    required=True, type=Path)
    p.add_argument("--tpm-matrix",       type=Path, default=None)
    p.add_argument("--te-bed",           type=Path, default=None)
    p.add_argument("--outdir",           required=True, type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    comp = load_bed(args.compartments_bed,
                    ["chrom", "start", "end", "label"])
    comp_by_chr = {c: g.reset_index(drop=True) for c, g in comp.groupby("chrom")}

    genes = load_bed(args.gene_bed,
                     ["chrom", "start", "end", "gene", "score", "strand"])
    sub_map = pd.read_csv(args.subgenome_map, sep="\t", header=None,
                          names=["gene", "subgenome"]).set_index("gene")

    genes["compartment"] = genes.apply(
        lambda r: assign_compartment(r, comp_by_chr), axis=1)
    genes = genes.merge(sub_map, left_on="gene", right_index=True, how="left")
    genes_out = args.outdir / "genes_by_compartment.tsv"
    genes.to_csv(genes_out, sep="\t", index=False)

    summary = (genes.groupby(["subgenome", "compartment"])
                    .size()
                    .reset_index(name="n_genes"))
    summary.to_csv(args.outdir / "compartment_gene_counts.tsv",
                   sep="\t", index=False)

    # TE coverage in A vs B
    if args.te_bed is not None:
        te = load_bed(args.te_bed, ["chrom", "start", "end", "name"])
        te["midpoint"] = (te["start"] + te["end"]) // 2
        te["compartment"] = te.apply(
            lambda r: assign_compartment(
                {"chrom": r["chrom"], "start": r["midpoint"],
                 "end": r["midpoint"] + 1}, comp_by_chr), axis=1)
        te_summary = te.groupby("compartment")["name"].count().reset_index(
            name="n_TEs")
        te_summary.to_csv(args.outdir / "compartment_TE_counts.tsv",
                          sep="\t", index=False)

    # per-tissue TPM by compartment (long format, easy to plot)
    if args.tpm_matrix is not None:
        tpm = pd.read_csv(args.tpm_matrix, sep="\t")
        tpm = tpm.rename(columns={tpm.columns[0]: "gene"})
        merged = genes[["gene", "subgenome", "compartment"]].merge(
            tpm, on="gene", how="inner")
        long = merged.melt(id_vars=["gene", "subgenome", "compartment"],
                           var_name="tissue", value_name="tpm")
        long.to_csv(args.outdir / "genes_compartment_tpm_long.tsv",
                    sep="\t", index=False)
        summary_tpm = (long.groupby(["subgenome", "compartment", "tissue"])
                            ["tpm"].mean().reset_index())
        summary_tpm.to_csv(args.outdir / "mean_tpm_by_subgenome_compartment.tsv",
                           sep="\t", index=False)

    with open(args.outdir / "ab_compartment_analysis.done", "w") as fh:
        json.dump({"n_genes": int(len(genes)),
                   "n_genes_with_compartment":
                       int((genes["compartment"] != "NA").sum())},
                  fh, indent=2)


if __name__ == "__main__":
    main()
