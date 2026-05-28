#!/usr/bin/env python3


from __future__ import annotations

import argparse
import json
from pathlib import Path

from scipy import stats


def load_counts(path: Path) -> tuple[int, int]:
    with open(path) as fh:
        d = json.load(fh)
    return int(d["n_D_silenced"]), int(d["n_S_silenced"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cave-json",    required=True, type=Path,
                   help="stats.json from cave species run")
    p.add_argument("--surface-json", required=True, type=Path,
                   help="stats.json from surface species run")
    p.add_argument("--out",          required=True, type=Path,
                   help="Output JSON with Fisher test results")
    args = p.parse_args()

    d_cave,    s_cave    = load_counts(args.cave_json)
    d_surface, s_surface = load_counts(args.surface_json)

    table = [[d_cave, s_cave],
             [d_surface, s_surface]]

    odds, p_value = stats.fisher_exact(table, alternative="two-sided")

    out = {
        "contingency_table": {
            "cave":    {"D_silenced": d_cave,    "S_silenced": s_cave},
            "surface": {"D_silenced": d_surface, "S_silenced": s_surface},
        },
        "fisher_odds_ratio":  float(odds),
        "fisher_p_twosided":  float(p_value),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"Fisher's exact test (cave vs surface): "
          f"OR = {odds:.3f}, P = {p_value:.2e}")


if __name__ == "__main__":
    main()
