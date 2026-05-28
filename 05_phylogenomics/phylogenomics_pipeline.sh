#!/usr/bin/env bash

set -euo pipefail

THREADS=16
while getopts "i:c:o:t:" opt; do
  case ${opt} in
    i) PROT_DIR=${OPTARG} ;;
    c) CDS_DIR=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${PROT_DIR:?Provide -i proteins_dir}"
: "${CDS_DIR:?Provide -c cds_dir}"
: "${OUTDIR:?Provide -o outdir}"

mkdir -p "$OUTDIR"

# -----------------------------------------------------------------------------
# 1. OrthoFinder
# -----------------------------------------------------------------------------
echo "[$(date +%T)] OrthoFinder ..."
orthofinder -f "$PROT_DIR" -t "$THREADS" -a "$THREADS" \
            -o "$OUTDIR/orthofinder"
OF_RES="$OUTDIR/orthofinder/Results_$(ls "$OUTDIR/orthofinder" | grep Results | head -1 | sed 's/Results_//')"
[[ -d "$OUTDIR/orthofinder" ]] && OF_RES=$(ls -d "$OUTDIR"/orthofinder/Results_* | head -1)
SC_DIR="$OF_RES/Single_Copy_Orthologue_Sequences"
echo "[INFO] single-copy OGs in: $SC_DIR"

# -----------------------------------------------------------------------------
# 2-3. MAFFT + pal2nal per OG
# -----------------------------------------------------------------------------
ALN_DIR="$OUTDIR/alignments"
CODON_DIR="$OUTDIR/codon_alignments"
mkdir -p "$ALN_DIR" "$CODON_DIR"

for og_fa in "$SC_DIR"/*.fa; do
  og=$(basename "$og_fa" .fa)
  mafft --auto --thread 4 "$og_fa" > "$ALN_DIR/${og}.aln" 2>/dev/null

  # Build matching CDS fasta from CDS_DIR
  python3 "$(dirname "$0")/extract_cds_for_og.py" \
    --aln-fa  "$og_fa" \
    --cds-dir "$CDS_DIR" \
    --out     "$ALN_DIR/${og}.cds.fa"

  pal2nal.pl "$ALN_DIR/${og}.aln" "$ALN_DIR/${og}.cds.fa" \
             -output paml -nogap > "$CODON_DIR/${og}.phy" \
             2> "$CODON_DIR/${og}.log" || true
done

# -----------------------------------------------------------------------------
# 4. Supermatrix concatenation (AMAS)
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Concatenating supermatrix ..."
AMAS.py concat \
  -i "$CODON_DIR"/*.phy \
  -f phylip -d dna \
  -t "$OUTDIR/supermatrix.phy" \
  -p "$OUTDIR/partitions.txt"

# -----------------------------------------------------------------------------
# 5. IQ-TREE ML inference
# -----------------------------------------------------------------------------
echo "[$(date +%T)] IQ-TREE ML ..."
iqtree2 \
  -s "$OUTDIR/supermatrix.phy" \
  -p "$OUTDIR/partitions.txt" \
  -m MFP \
  -B 1000 -alrt 1000 \
  -T "$THREADS" \
  --prefix "$OUTDIR/iqtree"

# -----------------------------------------------------------------------------
# 6. Divergence dating with treePL
# -----------------------------------------------------------------------------
echo "[$(date +%T)] treePL dating ..."

# Generate calibration block dynamically (user must edit taxon names below
# to match those in the IQ-TREE tree).
cat > "$OUTDIR/treepl.cfg" <<'CFG'
treefile = REPLACE_WITH_TREE
smooth   = 100
numsites = 0
mrca = NODE_DR_VS_CG Danio_rerio Cyprinus_carpio_D
min  = NODE_DR_VS_CG 32.81
max  = NODE_DR_VS_CG 36.84

mrca = NODE_DSUB_CARP_GOLDFISH Cyprinus_carpio_D Carassius_auratus_D
min  = NODE_DSUB_CARP_GOLDFISH 8.14
max  = NODE_DSUB_CARP_GOLDFISH 9.14

mrca = NODE_SSUB_CARP_GOLDFISH Cyprinus_carpio_S Carassius_auratus_S
min  = NODE_SSUB_CARP_GOLDFISH 7.74
max  = NODE_SSUB_CARP_GOLDFISH 8.69

outfile = treePL_dated.nwk
thorough
prime
CFG

sed -i "s|REPLACE_WITH_TREE|$OUTDIR/iqtree.treefile|" "$OUTDIR/treepl.cfg"
treePL "$OUTDIR/treepl.cfg" || echo "[WARN] treePL failed; edit calibrations in $OUTDIR/treepl.cfg and rerun."

echo "[$(date +%T)] Done."
