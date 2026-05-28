#!/usr/bin/env bash

set -euo pipefail

THREADS=32
while getopts "g:b:p:t:s:o:c:" opt; do
  case ${opt} in
    g) GENOME=${OPTARG} ;;
    b) RNABAM=${OPTARG} ;;
    p) PROTEINS=${OPTARG} ;;
    t) TRINITY=${OPTARG} ;;
    s) SPECIES=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    c) THREADS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${GENOME:?Provide -g <softmasked_genome.fa>}"
: "${RNABAM:?Provide -b <rnaseq.bam>}"
: "${PROTEINS:?Provide -p <proteins.fa>}"
: "${TRINITY:?Provide -t <trinity.fasta>}"
: "${SPECIES:?Provide -s <species_tag>}"
: "${OUTDIR:?Provide -o outdir}"

BRAKER_DIR="$OUTDIR/01_braker3"
PASA_DIR="$OUTDIR/02_pasa"
FUNC_DIR="$OUTDIR/03_functional"
mkdir -p "$BRAKER_DIR" "$PASA_DIR" "$FUNC_DIR"

# -----------------------------------------------------------------------------
# 1. BRAKER3: combines AUGUSTUS + GeneMark-ETP using RNA-seq + protein evidence;
#    BRAKER3 internally calls TSEBRA to combine the two prediction tracks.
# -----------------------------------------------------------------------------
echo "[$(date +%T)] BRAKER3..."
cd "$BRAKER_DIR"
braker.pl \
  --genome="$GENOME" \
  --species="$SPECIES" \
  --bam="$RNABAM" \
  --prot_seq="$PROTEINS" \
  --softmasking \
  --threads "$THREADS" \
  --gff3 \
  --workingdir="$BRAKER_DIR"

BRAKER_GFF="$BRAKER_DIR/braker.gff3"

# -----------------------------------------------------------------------------
# 2. PASA: align Trinity transcripts to the genome and update BRAKER3 models
#    with refined UTRs / alternative isoforms.
# -----------------------------------------------------------------------------
echo "[$(date +%T)] PASA..."
cd "$PASA_DIR"

# Step 2a: align + cluster transcripts
${PASAHOME:-/opt/PASApipeline}/Launch_PASA_pipeline.pl \
  -c "${PASAHOME:-/opt/PASApipeline}/pasa_conf/pasa.alignAssembly.Template.txt" \
  -C -R -g "$GENOME" -t "$TRINITY" \
  --ALIGNERS blat,gmap --CPU "$THREADS"

# Step 2b: load existing BRAKER3 annotation and update
${PASAHOME:-/opt/PASApipeline}/scripts/Load_Current_Gene_Annotations.dbi \
  -c "${PASAHOME:-/opt/PASApipeline}/pasa_conf/pasa.alignAssembly.Template.txt" \
  -g "$GENOME" -P "$BRAKER_GFF"

${PASAHOME:-/opt/PASApipeline}/Launch_PASA_pipeline.pl \
  -c "${PASAHOME:-/opt/PASApipeline}/pasa_conf/pasa.annotationCompare.Template.txt" \
  -A -g "$GENOME" -t "$TRINITY" --CPU "$THREADS"

FINAL_GFF="$PASA_DIR/pasa.assemblies.final.gff3"
FINAL_PEP="$PASA_DIR/pasa.assemblies.final.pep"

# Extract proteins (if not produced directly by PASA, regenerate with gffread)
gffread -y "$FINAL_PEP" -g "$GENOME" "$FINAL_GFF"

# -----------------------------------------------------------------------------
# 3. Functional annotation: InterProScan + eggNOG-mapper
# -----------------------------------------------------------------------------
echo "[$(date +%T)] InterProScan..."
cd "$FUNC_DIR"
# Strip stop codons (*) — InterProScan rejects them.
sed 's/\*$//' "$FINAL_PEP" > proteins.clean.fa

interproscan.sh \
  -i proteins.clean.fa \
  -f tsv,gff3,json \
  -goterms -pa -iprlookup \
  -cpu "$THREADS" \
  -b interproscan_out

echo "[$(date +%T)] eggNOG-mapper..."
emapper.py \
  -i proteins.clean.fa \
  --output eggnog_out \
  --output_dir "$FUNC_DIR" \
  --cpu "$THREADS" \
  -m diamond --itype proteins \
  --tax_scope Vertebrata \
  --go_evidence non-electronic

echo "[$(date +%T)] Gene annotation pipeline complete."
echo "Final GFF3:   $FINAL_GFF"
echo "Protein FA:   $FINAL_PEP"
echo "InterProScan: $FUNC_DIR/interproscan_out.tsv"
echo "eggNOG:       $FUNC_DIR/eggnog_out.emapper.annotations"
