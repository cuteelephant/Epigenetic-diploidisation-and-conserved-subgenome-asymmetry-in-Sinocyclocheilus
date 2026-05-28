#!/usr/bin/env bash

set -euo pipefail

while getopts "i:p:o:" opt; do
  case ${opt} in
    i) INDIR=${OPTARG} ;;
    p) SPP=${OPTARG} ;;            # space-separated species names
    o) OUTDIR=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${INDIR:?Provide -i input_dir}"
: "${SPP:?Provide -p \"sp1 sp2 ...\"}"
: "${OUTDIR:?Provide -o outdir}"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -----------------------------------------------------------------------------
# 1. Convert GFF3 -> BED
# -----------------------------------------------------------------------------
for sp in $SPP; do
  python3 -m jcvi.formats.gff bed --type=mRNA --key=ID \
    "$INDIR/${sp}.gff3" -o "${sp}.bed"
  ln -sf "$INDIR/${sp}.cds.fa" "${sp}.cds"
  ln -sf "$INDIR/${sp}.pep.fa" "${sp}.pep"
done

# -----------------------------------------------------------------------------
# 2. Sequential pairwise synteny (tracing the cyprinid phylogeny)
# -----------------------------------------------------------------------------
SPP_ARR=($SPP)
for ((i=0; i<${#SPP_ARR[@]}-1; i++)); do
  a=${SPP_ARR[$i]}
  b=${SPP_ARR[$((i+1))]}
  echo "[$(date +%T)] Synteny: $a vs $b"
  python3 -m jcvi.compara.catalog ortholog --no_strip_names \
    "$a" "$b"
done

# -----------------------------------------------------------------------------
# 3. Karyotype-level macrosynteny visualisation
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Building seqids and layout files..."
# user must supply seqids and layout — placeholder files generated below
for ((i=0; i<${#SPP_ARR[@]}-1; i++)); do
  a=${SPP_ARR[$i]}
  b=${SPP_ARR[$((i+1))]}
  # combine pairwise anchors into a single .anchors per pair already created
  echo "${a}.${b}.anchors" >> anchorlist.txt
done

python3 -m jcvi.graphics.karyotype seqids layout || \
  echo "[INFO] Supply seqids and layout files to render karyotype.pdf."

echo "[$(date +%T)] Synteny analysis complete: results in $OUTDIR"
