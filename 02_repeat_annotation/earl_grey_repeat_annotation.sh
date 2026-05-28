#!/usr/bin/env bash

set -euo pipefail

THREADS=32
while getopts "g:s:o:t:" opt; do
  case ${opt} in
    g) GENOME=${OPTARG} ;;
    s) SPECIES=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${GENOME:?Provide -g <genome.fa>}"
: "${SPECIES:?Provide -s <species_name>}"
: "${OUTDIR:?Provide -o outdir}"

mkdir -p "$OUTDIR"

earlGrey \
  -g "$GENOME" \
  -s "$SPECIES" \
  -o "$OUTDIR" \
  -r vertebrata \
  -t "$THREADS" \
  -e yes        # enable LTRfinder + RepeatCraft

echo "[$(date +%T)] Earl Grey outputs (BED/GFF + Kimura landscape) in: $OUTDIR/${SPECIES}_EarlGrey"
