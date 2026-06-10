#!/usr/bin/env bash
# Minimal usage example for euk-bin-refiner.
#
# This script is a template. Replace the paths with your own data before running.

set -euo pipefail

# === Paths ===
BINS_DIR="data/all_bins"                  # All your *.fa bin files together
BUSCO_CLASSIC="data/busco_classical"      # BUSCO results for CONCOCT/MaxBin/MetaBAT bins
BUSCO_REMAG="data/busco_remag"            # BUSCO results for REMAG bins
OUT_DIR="results/refined"

# === Parameters ===
ALPHA_D=5.0            # Duplication penalty
MIN_SCORE=0.05         # Bins below this can't win
MIN_CONTIGS=10         # Bins with fewer contigs are dropped

# === Run ===
python ../euk_bin_refiner.py \
    --bins-dir   "${BINS_DIR}" \
    --busco-dirs "${BUSCO_CLASSIC}" "${BUSCO_REMAG}" \
    --out-dir    "${OUT_DIR}" \
    --alpha-d    "${ALPHA_D}" \
    --min-score  "${MIN_SCORE}" \
    --min-contigs "${MIN_CONTIGS}" \
    --write-fasta

# === Inspect outputs ===
echo ""
echo "[INFO] Top winners (by BUSCO C):"
head -10 "${OUT_DIR}/winners.tsv"

echo ""
echo "[INFO] Cleaned FASTA bins:"
ls -1 "${OUT_DIR}/dereplicated_bins/" | head
