#!/usr/bin/env python3

from __future__ import annotations
import gzip
import io
import sys
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd


def smart_open(path: str | Path):
    """Open plain or gzipped text files transparently."""
    p = str(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
    return open(p, "r")


def load_tsv(path: str | Path, **kw) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kw)


def weighted_methylation(meth_c: np.ndarray, total_c: np.ndarray,
                         min_reads: int = 1) -> float:
    """
    Weighted methylation following Schultz et al. (2012):
        sum(meth_reads) / sum(total_reads)
    Returns NaN if no position has >= min_reads coverage.
    """
    mask = total_c >= min_reads
    if not mask.any() or total_c[mask].sum() == 0:
        return float("nan")
    return float(meth_c[mask].sum() / total_c[mask].sum())


def log2fc_tpm(tpm_a: float, tpm_b: float, pseudocount: float = 1.0) -> float:
    """Manuscript definition: log2((TPM_A + 1) / (TPM_B + 1))."""
    return float(np.log2((tpm_a + pseudocount) / (tpm_b + pseudocount)))


def fdr_bh(pvals: Iterable[float]) -> np.ndarray:
    """
    Benjamini–Hochberg adjusted p-values. NaN p-values are propagated.
    """
    p = np.asarray(list(pvals), dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = ~np.isnan(p)
    if not valid.any():
        return out
    pv = p[valid]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(1, n + 1))
    # Enforce monotonicity
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[valid] = res
    return out


def chunked(it: Iterable, n: int) -> Iterator[list]:
    buf = []
    for x in it:
        buf.append(x)
        if len(buf) == n:
            yield buf
            buf = []
    if buf:
        yield buf


def log(msg: str) -> None:
    print(f"[{Path(sys.argv[0]).name}] {msg}", file=sys.stderr)
