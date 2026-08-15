#!/usr/bin/env bash
# Generate N LOST synthetic centroid lists with known attitude.
# Usage: generate_lost_centroid_cases.sh N OUT_DIR DATABASE
set -euo pipefail

N="${1:?count}"
OUT_DIR="${2:?outdir}"
DATABASE="${3:?database}"
LOST_BIN="${4:?lost_binary}"
FOV="${5:-20}"

mkdir -p "${OUT_DIR}"
for i in $(seq 0 $((N - 1))); do
  SEED=$((20260802 + i))
  CASE=$(printf "lost_synth_%04d" "$i")
  TMP=$(mktemp)
  "${LOST_BIN}" pipeline \
    --generate 1 \
    --fov "${FOV}" \
    --generate-seed "${SEED}" \
    --generate-random-attitudes true \
    --generate-centroids-only true \
    --generate-spread-stddev 1 \
    --generate-read-noise-stddev 0.05 \
    --database "${DATABASE}" \
    --star-id-algo py \
    --angular-tolerance 0.05 \
    --false-stars-estimate 1000 \
    --max-mismatch-probability 0.0001 \
    --attitude-algo dqm \
    --print-expected-attitude - \
    --print-input-centroids "${TMP}" \
    > "${OUT_DIR}/${CASE}.attitude.txt" 2>/dev/null || true

  # Convert LOST print-input-centroids to TSV: x y intensity
  python3 - "${TMP}" "${OUT_DIR}/${CASE}.tsv" "${OUT_DIR}/${CASE}.json" "${OUT_DIR}/${CASE}.attitude.txt" <<'PY'
import json, re, sys
from pathlib import Path
src, tsv, js, att = map(Path, sys.argv[1:])
text = src.read_text(encoding="utf-8", errors="replace")
xs, ys = {}, {}
for line in text.splitlines():
    m = re.match(r"input_centroid_(\d+)_x\s+(\S+)", line.strip())
    if m:
        xs[int(m.group(1))] = float(m.group(2))
        continue
    m = re.match(r"input_centroid_(\d+)_y\s+(\S+)", line.strip())
    if m:
        ys[int(m.group(1))] = float(m.group(2))
att_text = att.read_text(encoding="utf-8", errors="replace")
truth = {}
for key, dest in (("expected_attitude_ra", "ra_deg"), ("expected_attitude_de", "dec_deg"), ("expected_attitude_roll", "roll_deg")):
    m = re.search(rf"{key}\s+(\S+)", att_text)
    if m:
        truth[dest] = float(m.group(1))
indices = sorted(set(xs) & set(ys))
lines = ["# x y intensity"]
centroids = []
for rank, idx in enumerate(indices):
    intensity = 1000 - rank
    lines.append(f"{xs[idx]:.6f} {ys[idx]:.6f} {intensity}")
    centroids.append({"x": xs[idx], "y": ys[idx], "intensity": intensity, "source": "lost_generate"})
tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
payload = {
    "case_id": tsv.stem,
    "provenance": "LOST --generate --generate-centroids-only",
    "truth": truth,
    "image_shape": [1024, 1024],
    "extractors": {"lost_generate": centroids},
}
js.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(f"wrote {tsv.name} n={len(centroids)} truth={truth}")
PY
  rm -f "${TMP}"
done
