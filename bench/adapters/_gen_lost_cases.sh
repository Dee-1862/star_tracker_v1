#!/usr/bin/env bash
set -euo pipefail

export LOST_BSC_PATH=/mnt/c/Users/Reddy/Desktop/st/star_tracker_suite/bench/data/hipparcos_for_lost.tsv
LOST=/mnt/c/Users/Reddy/Desktop/st/star_tracker_suite/benchmarking/lost_baseline/lost
DB=/mnt/c/Users/Reddy/Desktop/st/star_tracker_suite/bench/data/lost_database_hip.dat
OUT=/mnt/c/Users/Reddy/Desktop/st/star_tracker_suite/bench/centroids/lost_native
mkdir -p "${OUT}"

for i in 0 1 2 3 4; do
  SEED=$((20260802 + i))
  CASE=$(printf 'lost_synth_%04d' "$i")
  ATT="${OUT}/${CASE}.attitude.txt"
  CENT="${OUT}/${CASE}.centroids.txt"
  RA=$((20 + i * 30))
  DE=$((i * 10 - 20))
  ROLL=$((i * 40))
  "${LOST}" pipeline \
    --generate 1 \
    --fov 20 \
    --generate-seed "${SEED}" \
    --generate-ra "${RA}" \
    --generate-de "${DE}" \
    --generate-roll "${ROLL}" \
    --generate-centroids-only true \
    --generate-spread-stddev 1 \
    --database "${DB}" \
    --star-id-algo py \
    --attitude-algo dqm \
    --print-expected-attitude - \
    --print-input-centroids "${CENT}" \
    > "${ATT}" 2>/dev/null || true
  echo "CASE ${CASE}"
  grep expected_attitude_ra "${ATT}" || true
  grep num_input_centroids "${CENT}" | head -1 || true
done
