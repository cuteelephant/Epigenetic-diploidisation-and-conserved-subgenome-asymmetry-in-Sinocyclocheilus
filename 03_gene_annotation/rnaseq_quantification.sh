#!/usr/bin/env bash

set -euo pipefail

THREADS=16
while getopts "g:a:1:2:s:o:t:" opt; do
  case ${opt} in
    g) GENOME=${OPTARG} ;;
    a) GFF=${OPTARG} ;;
    1) R1=${OPTARG} ;;
    2) R2=${OPTARG} ;;
    s) SAMPLE=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${GENOME:?Provide -g <genome.fa>}"
: "${GFF:?Provide -a <annotation.gff3>}"
: "${R1:?Provide -1 R1}"; : "${R2:?Provide -2 R2}"
: "${SAMPLE:?Provide -s sample_tag}"
: "${OUTDIR:?Provide -o outdir}"

mkdir -p "$OUTDIR"
INDEX="$OUTDIR/hisat2_index"
BAM="$OUTDIR/${SAMPLE}.sorted.bam"

# Build HISAT2 index once (skip if present)
if [ ! -f "${INDEX}.1.ht2" ]; then
  echo "[$(date +%T)] Building HISAT2 index..."
  hisat2-build -p "$THREADS" "$GENOME" "$INDEX"
fi

echo "[$(date +%T)] HISAT2 align: $SAMPLE"
hisat2 -x "$INDEX" -1 "$R1" -2 "$R2" \
       --dta -p "$THREADS" --no-unal \
  | samtools sort -@ "$THREADS" -o "$BAM" -
samtools index "$BAM"

echo "[$(date +%T)] StringTie quantification (-eB)..."
mkdir -p "$OUTDIR/stringtie_${SAMPLE}"
stringtie "$BAM" \
  -G "$GFF" -eB \
  -o "$OUTDIR/stringtie_${SAMPLE}/${SAMPLE}.gtf" \
  -A "$OUTDIR/stringtie_${SAMPLE}/${SAMPLE}.gene_abund.tab" \
  -p "$THREADS"

echo "[$(date +%T)] Done. Gene-level TPM: $OUTDIR/stringtie_${SAMPLE}/${SAMPLE}.gene_abund.tab"
