#!/usr/bin/env bash

set -euo pipefail

THREADS=16
SAMPLE="sample"

while getopts "g:1:2:o:t:s:" opt; do
  case ${opt} in
    g) GENOME=${OPTARG} ;;
    1) FQ1=${OPTARG} ;;
    2) FQ2=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    s) SAMPLE=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${GENOME:?Provide -g <ref.fa>}"
: "${FQ1:?Provide -1 R1.fq.gz}"
: "${FQ2:?Provide -2 R2.fq.gz}"
: "${OUTDIR:?Provide -o outdir}"

MAP_DIR="${OUTDIR}/01_align/${SAMPLE}"
MAT_DIR="${OUTDIR}/02_matrix/${SAMPLE}"
COMP_DIR="${OUTDIR}/03_compartments/${SAMPLE}"
TAD_DIR="${OUTDIR}/04_tads/${SAMPLE}"
mkdir -p "$MAP_DIR" "$MAT_DIR" "$COMP_DIR" "$TAD_DIR"

# -----------------------------------------------------------------------------
# 1. Index genome (only once)
# -----------------------------------------------------------------------------
if [[ ! -f "${GENOME}.bwt" ]]; then
  bwa index "$GENOME"
fi

# -----------------------------------------------------------------------------
# 2. Hi-C alignment with the Methods-mandated parameters
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Aligning Hi-C R1 ..."
bwa mem -A1 -B4 -E50 -L0 -t "$THREADS" "$GENOME" "$FQ1" \
  | samtools view -bS - > "$MAP_DIR/${SAMPLE}_R1.bam"

echo "[$(date +%T)] Aligning Hi-C R2 ..."
bwa mem -A1 -B4 -E50 -L0 -t "$THREADS" "$GENOME" "$FQ2" \
  | samtools view -bS - > "$MAP_DIR/${SAMPLE}_R2.bam"

# -----------------------------------------------------------------------------
# 3. Build matrices with HiCExplorer at 100kb and 40kb resolutions
# -----------------------------------------------------------------------------
echo "[$(date +%T)] hicBuildMatrix at 100kb ..."
hicBuildMatrix \
  --samFiles "$MAP_DIR/${SAMPLE}_R1.bam" "$MAP_DIR/${SAMPLE}_R2.bam" \
  --binSize  100000 \
  --restrictionSequence GATC \
  --danglingSequence    GATC \
  --threads  "$THREADS" \
  --outBam   "$MAT_DIR/${SAMPLE}_100kb.bam" \
  --outFileName "$MAT_DIR/${SAMPLE}_100kb.h5" \
  --QCfolder "$MAT_DIR/${SAMPLE}_100kb_QC"

echo "[$(date +%T)] hicBuildMatrix at 40kb ..."
hicBuildMatrix \
  --samFiles "$MAP_DIR/${SAMPLE}_R1.bam" "$MAP_DIR/${SAMPLE}_R2.bam" \
  --binSize  40000 \
  --restrictionSequence GATC \
  --danglingSequence    GATC \
  --threads  "$THREADS" \
  --outBam   "$MAT_DIR/${SAMPLE}_40kb.bam" \
  --outFileName "$MAT_DIR/${SAMPLE}_40kb.h5" \
  --QCfolder "$MAT_DIR/${SAMPLE}_40kb_QC"

# -----------------------------------------------------------------------------
# 4. ICE correction
# -----------------------------------------------------------------------------
echo "[$(date +%T)] ICE-correcting matrices ..."
hicCorrectMatrix correct \
  --matrix       "$MAT_DIR/${SAMPLE}_100kb.h5" \
  --filterThreshold -1.5 5 \
  --outFileName  "$MAT_DIR/${SAMPLE}_100kb_corr.h5"

hicCorrectMatrix correct \
  --matrix       "$MAT_DIR/${SAMPLE}_40kb.h5" \
  --filterThreshold -1.5 5 \
  --outFileName  "$MAT_DIR/${SAMPLE}_40kb_corr.h5"

# -----------------------------------------------------------------------------
# 5. A/B compartments via CALDER (requires .cool input)
#    Conversion h5 -> cool first.
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Converting to .cool for CALDER ..."
hicConvertFormat --matrices "$MAT_DIR/${SAMPLE}_100kb_corr.h5" \
                 --outFileName "$MAT_DIR/${SAMPLE}_100kb_corr.cool" \
                 --inputFormat h5 --outputFormat cool

# CALDER call (R package). Wrapper Rscript is provided alongside:
echo "[$(date +%T)] Calling A/B compartments via CALDER ..."
Rscript --vanilla "$(dirname "$0")/calder_compartments.R" \
        "$MAT_DIR/${SAMPLE}_100kb_corr.cool" \
        100000 \
        "$COMP_DIR"

# -----------------------------------------------------------------------------
# 6. TAD calling
# -----------------------------------------------------------------------------
echo "[$(date +%T)] hicFindTADs ..."
hicFindTADs \
  --matrix      "$MAT_DIR/${SAMPLE}_40kb_corr.h5" \
  --outPrefix   "$TAD_DIR/${SAMPLE}_TAD" \
  --minDepth    120000 \
  --maxDepth    480000 \
  --step         40000 \
  --thresholdComparisons 0.01 \
  --delta       0.01 \
  --correctForMultipleTesting fdr \
  --numberOfProcessors "$THREADS"

echo "[$(date +%T)] Done. A/B in: $COMP_DIR ; TADs in: $TAD_DIR"
