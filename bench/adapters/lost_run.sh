#!/usr/bin/env bash
# Run LOST full pipeline on a grayscale PNG. Requires LOST built on Linux/WSL.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LOST_BIN="${ROOT}/benchmarking/lost_baseline/lost"
DATABASE="${ROOT}/bench/data/lost_database_20deg.dat"
PNG="${1:?usage: lost_run.sh image.png [fov_deg]}"
FOV="${2:-20}"

if [[ ! -x "${LOST_BIN}" ]]; then
  echo "LOST binary not found. Build with: cd benchmarking/lost_baseline && make LOST_DISABLE_ASAN=1 CXX=g++-12" >&2
  exit 1
fi

if [[ ! -f "${DATABASE}" ]]; then
  echo "LOST database not found at ${DATABASE}. See benchmarking/lost_baseline/BASELINE.md" >&2
  exit 1
fi

exec "${LOST_BIN}" pipeline \
  --png "${PNG}" \
  --fov "${FOV}" \
  --centroid-algo cog \
  --centroid-mag-filter 5 \
  --database "${DATABASE}" \
  --star-id-algo py \
  --angular-tolerance 0.05 \
  --false-stars-estimate 1000 \
  --max-mismatch-probability 0.0001 \
  --attitude-algo dqm \
  --print-attitude - \
  --print-speed -
