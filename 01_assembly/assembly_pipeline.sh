#!/usr/bin/env bash

set -euo pipefail

THREADS=32
N_CHRS=50
while getopts "p:1:2:o:t:n:" opt; do
  case ${opt} in
    p) HIFI=${OPTARG} ;;
    1) HIC1=${OPTARG} ;;
    2) HIC2=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    n) N_CHRS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

: "${HIFI:?Provide -p <pacbio_hifi.fq.gz>}"
: "${HIC1:?Provide -1 <hic_R1.fq.gz>}"
: "${HIC2:?Provide -2 <hic_R2.fq.gz>}"
: "${OUTDIR:?Provide -o outdir}"

HIFIASM_DIR="$OUTDIR/01_hifiasm"
PURGE_DIR="$OUTDIR/02_purge_dups"
HAPHIC_DIR="$OUTDIR/03_haphic"
BUSCO_DIR="$OUTDIR/04_busco"
mkdir -p "$HIFIASM_DIR" "$PURGE_DIR" "$HAPHIC_DIR" "$BUSCO_DIR"

# -----------------------------------------------------------------------------
# 1. Hifiasm
# -----------------------------------------------------------------------------
echo "[$(date +%T)] Hifiasm..."
hifiasm -o "$HIFIASM_DIR/asm" -t "$THREADS" "$HIFI"
awk '/^S/{print ">"$2"\n"$3}' "$HIFIASM_DIR/asm.bp.p_ctg.gfa" \
  > "$HIFIASM_DIR/asm.p_ctg.fa"

# -----------------------------------------------------------------------------
# 2. purge_dups
# -----------------------------------------------------------------------------
echo "[$(date +%T)] purge_dups..."
cd "$PURGE_DIR"
minimap2 -t "$THREADS" -xmap-pb "$HIFIASM_DIR/asm.p_ctg.fa" "$HIFI" \
  | gzip -c > reads.paf.gz
pbcstat reads.paf.gz
calcuts PB.stat > cutoffs 2> calcuts.log
split_fa "$HIFIASM_DIR/asm.p_ctg.fa" > asm.split
minimap2 -t "$THREADS" -xasm5 -DP asm.split asm.split \
  | gzip -c > asm.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov asm.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed "$HIFIASM_DIR/asm.p_ctg.fa"
mv purged.fa "$PURGE_DIR/asm.purged.fa"

# -----------------------------------------------------------------------------
# 3. HapHiC scaffolding (Hi-C)
# -----------------------------------------------------------------------------
echo "[$(date +%T)] HapHiC..."
cd "$HAPHIC_DIR"
ln -sf "$PURGE_DIR/asm.purged.fa" .

# index for BWA-MEM
bwa index asm.purged.fa
bwa mem -5SP -t "$THREADS" asm.purged.fa "$HIC1" "$HIC2" \
  | samblaster \
  | samtools view -F 3340 -hb - \
  > hic.bam

haphic pipeline asm.purged.fa hic.bam "$N_CHRS" --threads "$THREADS" \
  --RE GATC --remove_allelic_links 2

# Convert to .hic for Juicebox manual review
haphic juicer scaffolds.fa hic.bam

# -----------------------------------------------------------------------------
# 4. BUSCO
# -----------------------------------------------------------------------------
echo "[$(date +%T)] BUSCO..."
cd "$BUSCO_DIR"
busco -i "$HAPHIC_DIR/scaffolds.fa" \
      -l actinopterygii_odb10 \
      -m genome --cpu "$THREADS" -o busco_out

echo "[$(date +%T)] Pipeline complete."
echo "Final scaffolds: $HAPHIC_DIR/scaffolds.fa"
