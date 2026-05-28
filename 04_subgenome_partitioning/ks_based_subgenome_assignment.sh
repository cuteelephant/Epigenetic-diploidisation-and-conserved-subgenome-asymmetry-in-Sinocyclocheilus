#!/usr/bin/env bash

set -euo pipefail

THREADS=8
while getopts "q:t:G:T:o:p:" opt; do
  case ${opt} in
    q) QPROT=${OPTARG} ;;
    t) TPROT=${OPTARG} ;;
    G) QGFF=${OPTARG} ;;
    T) TGFF=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
    p) THREADS=${OPTARG} ;;
    *) echo "Unknown option"; exit 1 ;;
  esac
done

mkdir -p "$OUTDIR"

# -----------------------------------------------------------------------------
# 1. Protein alignment with DIAMOND
# -----------------------------------------------------------------------------
diamond makedb --in "$TPROT" -d "$OUTDIR/puntius"
diamond blastp \
  --query "$QPROT" \
  --db    "$OUTDIR/puntius" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --evalue 1e-5 \
  --threads "$THREADS" \
  --out "$OUTDIR/diamond_q_vs_puntius.tsv"

# -----------------------------------------------------------------------------
# 2. WGDI: collinearity blocks, Ks and integration
# -----------------------------------------------------------------------------
# WGDI uses a config file; we generate a minimal one inline.
CONF="$OUTDIR/wgdi.conf"
cat > "$CONF" <<EOF
[ks]
cds_file = ${QPROT%.faa}.cds
pep_file = ${QPROT}
align_software = mafft
pairs_file = ${OUTDIR}/diamond_q_vs_puntius.tsv
ks_file = ${OUTDIR}/q_vs_puntius.ks

[collinearity]
gff1 = ${QGFF}
gff2 = ${TGFF}
lens1 = ${OUTDIR}/q.lens
lens2 = ${OUTDIR}/t.lens
blast = ${OUTDIR}/diamond_q_vs_puntius.tsv
collinearity = ${OUTDIR}/q_vs_puntius.collinearity
score = 100
evalue = 1e-5
grading = 50,40,25
mg = 25,25
pvalue = 1
EOF

wgdi -icl "$CONF"
wgdi -ks  "$CONF"
wgdi -bi  "$CONF"
wgdi -c   "$CONF"     # 1:1 syntenic filtering
wgdi -bk  "$CONF"     # Ks-coloured dotplot, used for validation

# -----------------------------------------------------------------------------
# 3. Per-chromosome median Ks summary  -> input to assign SubD/SubS
# -----------------------------------------------------------------------------
python3 - <<'PY'
import os, sys, pandas as pd

outdir = os.environ.get("OUTDIR", ".")
ks_path = os.path.join(outdir, "q_vs_puntius.ks")
if not os.path.exists(ks_path):
    sys.exit("[ERROR] missing Ks output: " + ks_path)

ks = pd.read_csv(ks_path, sep="\t")
# WGDI Ks columns: id1, id2, ka_NG86, ks_NG86, ka_YN00, ks_YN00 (may vary)
chr_col_candidates = [c for c in ks.columns if "chr" in c.lower() or c == "id1"]
chr_col = chr_col_candidates[0] if chr_col_candidates else "id1"
ks_col = "ks_YN00" if "ks_YN00" in ks.columns else ks.columns[3]

ks["chrom_q"] = ks[chr_col].str.split(".").str[0]
med = (ks.groupby("chrom_q")[ks_col]
         .median()
         .reset_index()
         .rename(columns={ks_col: "median_Ks"}))
med = med.sort_values("median_Ks")
# split at median
cutoff = med["median_Ks"].median()
med["assignment"] = ["SubD" if k <= cutoff else "SubS"
                     for k in med["median_Ks"]]
med.to_csv(os.path.join(outdir, "chromosome_Ks_assignment.tsv"),
           sep="\t", index=False)
print(med.to_string(index=False))
PY

echo "[$(date +%T)] Ks-based subgenome assignment written to $OUTDIR/chromosome_Ks_assignment.tsv"
