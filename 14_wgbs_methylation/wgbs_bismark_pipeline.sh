#!/usr/bin/env bash

set -euo pipefail

THREADS=8
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

REF_DIR="$(dirname "$GENOME")"
TRIM_DIR="${OUTDIR}/01_trim_galore/${SAMPLE}"
MAP_DIR="${OUTDIR}/02_bismark/${SAMPLE}"
EXTRACT_DIR="${OUTDIR}/03_methylation/${SAMPLE}"
mkdir -p "$TRIM_DIR" "$MAP_DIR" "$EXTRACT_DIR"

# -----------------------------------------------------------------------------
# 1. Genome preparation (only once per reference)
# -----------------------------------------------------------------------------
if [[ ! -d "${REF_DIR}/Bisulfite_Genome" ]]; then
  echo "[$(date +%T)] Preparing bisulfite-converted genome ..."
  bismark_genome_preparation --bowtie2 "$REF_DIR"
fi

# -----------------------------------------------------------------------------
# 2. Adapter / quality trimming
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Trim Galore on ${SAMPLE} ..."
trim_galore --paired --cores "$THREADS" \
            --output_dir "$TRIM_DIR" \
            "$FQ1" "$FQ2"

R1_TRIM=$(ls "$TRIM_DIR"/*val_1.fq.gz | head -1)
R2_TRIM=$(ls "$TRIM_DIR"/*val_2.fq.gz | head -1)

# -----------------------------------------------------------------------------
# 3. Bismark alignment (bowtie2 backend)
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Bismark mapping ${SAMPLE} ..."
bismark --bowtie2 --multicore "$((THREADS / 2))" \
        --genome "$REF_DIR" \
        -1 "$R1_TRIM" -2 "$R2_TRIM" \
        -o "$MAP_DIR"

BAM=$(ls "$MAP_DIR"/*_pe.bam | head -1)

# -----------------------------------------------------------------------------
# 4. Deduplication
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Deduplicating ..."
deduplicate_bismark --paired --bam --output_dir "$MAP_DIR" "$BAM"
DEDUP_BAM=$(ls "$MAP_DIR"/*deduplicated.bam | head -1)

# -----------------------------------------------------------------------------
# 5. Methylation extraction (CpG / CHG / CHH)
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Extracting methylation ..."
bismark_methylation_extractor --paired-end --comprehensive --bedGraph \
                              --cytosine_report --genome_folder "$REF_DIR" \
                              --multicore "$((THREADS / 2))" \
                              --output "$EXTRACT_DIR" \
                              "$DEDUP_BAM"

echo "[$(date +%T)] Done. Outputs in: $EXTRACT_DIR"
