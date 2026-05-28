#!/usr/bin/env bash

set -euo pipefail

while getopts "g:c:o:" opt; do
  case ${opt} in
    g) GENOME=${OPTARG} ;;
    c) PAIRS=${OPTARG}  ;;
    o) OUTDIR=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${GENOME:?Provide -g <genome.fa>}"
: "${PAIRS:?Provide -c <chromosome_pairs.tsv>}"
: "${OUTDIR:?Provide -o outdir}"

mkdir -p "$OUTDIR"

subphaser \
  -i "$GENOME" \
  -c "$PAIRS" \
  -k 15 \
  -q 200 \
  -f 2 \
  -bg_bin 100000 \
  -tmp_dir "$OUTDIR/tmp" \
  -o "$OUTDIR" \
  -tre

echo "[$(date +%T)] SubPhaser outputs in $OUTDIR"
echo "  -- chromosome_assignments.tsv : SubD/SubS labels"
echo "  -- pca.pdf, dendrogram.pdf    : validation panels for Fig. 1c,d"
