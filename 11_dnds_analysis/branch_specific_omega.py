#!/usr/bin/env python3


from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CODEML_CTL_TEMPLATE = """\
      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 2
      NSsites = 0
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 8
       getSE  = 0
 RateAncestor = 0
    cleandata = 1
"""

OMEGA_LINE_RE = re.compile(r"w \(dN/dS\) for branches:\s+([0-9eE\.\-\+ ]+)")


# ---------------------------------------------------------------------------
# Tree labelling
# ---------------------------------------------------------------------------
def label_tree_for_og(tree_str: str, subd_tip: str, subs_tip: str) -> str:
    """Insert '#1' on the branch leading to subd_tip and '#2' on subs_tip.

    Uses a simple text substitution that matches the tip label followed by
    ':' (branch length). codeml's tree-label syntax accepts this on terminal
    branches; for internal branches the user should supply a tree where the
    intended ancestor is already represented as a labelled internal node.

    NOTE: For internal-branch labelling you may want to use ete3 / dendropy.
    Here we keep dependencies minimal.
    """
    out = re.sub(rf"\b{re.escape(subd_tip)}\b(?=:)", f"{subd_tip} #1", tree_str)
    out = re.sub(rf"\b{re.escape(subs_tip)}\b(?=:)", f"{subs_tip} #2", out)
    return out


# ---------------------------------------------------------------------------
# Per-OG codeml driver
# ---------------------------------------------------------------------------
def run_one_og(args_tuple) -> dict:
    aln_path, tree_str, codeml_bin, subd_tip, subs_tip, workdir = args_tuple
    og_id = aln_path.stem
    og_work = workdir / og_id
    og_work.mkdir(parents=True, exist_ok=True)

    seqfile = og_work / "aln.phy"
    shutil.copy(aln_path, seqfile)

    treefile = og_work / "tree.nwk"
    with open(treefile, "w") as fh:
        fh.write(label_tree_for_og(tree_str, subd_tip, subs_tip))

    outfile = og_work / "mlc.txt"
    ctl_file = og_work / "codeml.ctl"
    with open(ctl_file, "w") as fh:
        fh.write(CODEML_CTL_TEMPLATE.format(
            seqfile=seqfile.name,
            treefile=treefile.name,
            outfile=outfile.name))

    try:
        subprocess.run([codeml_bin, ctl_file.name],
                       cwd=og_work, check=True,
                       capture_output=True, timeout=3600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        return {"OG": og_id, "status": f"codeml_failed:{type(e).__name__}",
                "omega_D": np.nan, "omega_S": np.nan}

    # parse omega values
    try:
        text = outfile.read_text()
    except FileNotFoundError:
        return {"OG": og_id, "status": "no_outfile",
                "omega_D": np.nan, "omega_S": np.nan}

    m = OMEGA_LINE_RE.search(text)
    if not m:
        return {"OG": og_id, "status": "parse_failed",
                "omega_D": np.nan, "omega_S": np.nan}
    parts = [float(x) for x in m.group(1).split() if x]
    # codeml reports omega for background then foregrounds in label order
    # we read the last two as foregrounds (#1, #2)
    if len(parts) < 3:
        return {"OG": og_id, "status": "too_few_omegas",
                "omega_D": np.nan, "omega_S": np.nan}
    omega_d, omega_s = parts[-2], parts[-1]
    return {"OG": og_id, "status": "ok",
            "omega_D": omega_d, "omega_S": omega_s}


# ---------------------------------------------------------------------------
# Filtering and tests
# ---------------------------------------------------------------------------
def filter_omegas(df: pd.DataFrame,
                  low: float = 5e-4, high: float = 10.0) -> pd.DataFrame:
    return df[(df["omega_D"] > low) & (df["omega_D"] < high) &
              (df["omega_S"] > low) & (df["omega_S"] < high) &
              df["omega_D"].notna() & df["omega_S"].notna()].copy()


def selection_regime(omega: float) -> str:
    if omega < 0.1:
        return "strong"
    if omega < 0.5:
        return "moderate"
    if omega < 1.0:
        return "relaxed"
    return "positive"


def run_tests(df: pd.DataFrame) -> dict:
    if df.empty:
        return {"n": 0}

    paired_two = stats.wilcoxon(df["omega_D"], df["omega_S"],
                                alternative="two-sided")
    paired_one = stats.wilcoxon(df["omega_D"], df["omega_S"],
                                alternative="greater")

    df["regime_D"] = df["omega_D"].map(selection_regime)
    df["regime_S"] = df["omega_S"].map(selection_regime)
    regimes = ["strong", "moderate", "relaxed", "positive"]
    ct = np.array([[(df["regime_D"] == r).sum() for r in regimes],
                   [(df["regime_S"] == r).sum() for r in regimes]])
    chi2, p_chi, dof, _ = stats.chi2_contingency(ct)

    # McNemar on omega>1 status
    a = ((df["omega_D"] >  1) & (df["omega_S"] <= 1)).sum()
    b = ((df["omega_D"] <= 1) & (df["omega_S"] >  1)).sum()
    # exact McNemar
    if a + b > 0:
        mc_p = stats.binomtest(int(min(a, b)), int(a + b), p=0.5,
                               alternative="two-sided").pvalue
    else:
        mc_p = float("nan")

    return {
        "n_passed":            int(len(df)),
        "median_omega_D":      float(df["omega_D"].median()),
        "median_omega_S":      float(df["omega_S"].median()),
        "paired_wilcoxon_two_sided": {
            "W": float(paired_two.statistic),
            "p": float(paired_two.pvalue),
        },
        "paired_wilcoxon_D_gt_S": {
            "W": float(paired_one.statistic),
            "p": float(paired_one.pvalue),
        },
        "regime_counts": {
            "SubD": dict(zip(regimes, ct[0].tolist())),
            "SubS": dict(zip(regimes, ct[1].tolist())),
        },
        "chi2": {
            "statistic": float(chi2),
            "df":        int(dof),
            "p":         float(p_chi),
        },
        "mcnemar_omega_gt_1": {
            "a_only_D":  int(a),
            "b_only_S":  int(b),
            "p":         float(mc_p) if not np.isnan(mc_p) else None,
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--aln-dir",  required=True, type=Path)
    p.add_argument("--tree",     required=True, type=Path,
                   help="Newick tree shared across OGs")
    p.add_argument("--codeml",   default="codeml",
                   help="codeml executable (PAML)")
    p.add_argument("--subd-tip", required=True,
                   help="tip label whose ancestral edge is SubD (#1)")
    p.add_argument("--subs-tip", required=True,
                   help="tip label whose ancestral edge is SubS (#2)")
    p.add_argument("--workdir",  required=True, type=Path)
    p.add_argument("--threads",  type=int, default=4)
    p.add_argument("--outdir",   required=True, type=Path)
    args = p.parse_args()

    args.workdir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    tree_str = args.tree.read_text().strip()
    aln_files = sorted(args.aln_dir.glob("*.phy")) + \
                sorted(args.aln_dir.glob("*.fa")) + \
                sorted(args.aln_dir.glob("*.fasta"))
    if not aln_files:
        sys.exit(f"[ERROR] no alignments found in {args.aln_dir}")

    tasks = [(f, tree_str, args.codeml, args.subd_tip, args.subs_tip,
              args.workdir) for f in aln_files]
    with mp.Pool(args.threads) as pool:
        records = pool.map(run_one_og, tasks)

    table = pd.DataFrame(records)
    table.to_csv(args.outdir / "omega_table.tsv", sep="\t", index=False)

    passed = filter_omegas(table)
    passed.to_csv(args.outdir / "omega_passed.tsv", sep="\t", index=False)

    stats_dict = run_tests(passed)
    with open(args.outdir / "omega_stats.json", "w") as fh:
        json.dump(stats_dict, fh, indent=2)
    print(json.dumps(stats_dict, indent=2))


if __name__ == "__main__":
    main()
