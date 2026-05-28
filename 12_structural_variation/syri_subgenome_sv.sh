#!/usr/bin/env bash

set -euo pipefail

THREADS=8
while getopts "d:s:o:t:" opt; do
  case ${opt} in
    d) SUBD=${OPTARG} ;;
    s) SUBS=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${SUBD:?Provide -d <subD.fa>}"
: "${SUBS:?Provide -s <subS.fa>}"
: "${OUTDIR:?Provide -o outdir}"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -----------------------------------------------------------------------------
# 1. Whole-genome alignment with MUMmer4
# -----------------------------------------------------------------------------
echo "[$(date +%T)] nucmer ..."
nucmer --maxmatch -c 100 -b 500 -l 50 -t "$THREADS" \
       "$SUBD" "$SUBS" -p subD_vs_subS

echo "[$(date +%T)] delta-filter ..."
delta-filter -m -i 90 -l 100 subD_vs_subS.delta > subD_vs_subS.filt.delta
show-coords -THrd subD_vs_subS.filt.delta > subD_vs_subS.coords

# -----------------------------------------------------------------------------
# 2. SyRI structural-variation calling
# -----------------------------------------------------------------------------
echo "[$(date +%T)] SyRI ..."
syri \
  -c subD_vs_subS.coords \
  -d subD_vs_subS.filt.delta \
  -r "$SUBD" -q "$SUBS" \
  --nc "$THREADS" \
  --prefix subD_vs_subS_

# -----------------------------------------------------------------------------
# 3. Plot
# -----------------------------------------------------------------------------
echo "[$(date +%T)] plotsr ..."
plotsr \
  --sr subD_vs_subS_syri.out \
  --genomes <(printf "%s\t%s\n%s\t%s\n" SubD "$SUBD" SubS "$SUBS") \
  -o subD_vs_subS_plotsr.pdf || \
  echo "[INFO] plotsr requires a genomes.txt file; sample inline failed but main SyRI output is fine."

echo "[$(date +%T)] SV results: $OUTDIR/subD_vs_subS_syri.out"
