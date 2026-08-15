# Three-way star tracker benchmark harness

Compares three solvers on the same inputs:

| Solver | Location | Algorithm | Role |
|--------|----------|-----------|------|
| **tetra3** | `benchmarking/tetra3_baseline/` | 4-star edge-ratio hash + binomial verify | Reference (ESA) |
| **LOST** | `benchmarking/lost_baseline/` | Pyramid + k-vector + DQM attitude | Flight-realistic C++ baseline |
| **Ours** | `flight_software/` | Pair voting (today); pyramid + temporal LIS (planned) | Candidate |

## Design principle: decouple centroiding from star-ID

A solver can win because its centroiding is better, not because its pattern matching is better.
The harness therefore supports two modes:

1. **Full pipeline** — image in, attitude out (current `run_comparison.py` + LOST wrapper).
2. **Star-ID only** — frozen centroid JSON in, attitude/IDs out (`adapters/*_from_centroids`).

Frozen centroids live in `centroids/`. Export them once with `export_centroids.py`.

## Directory layout

```
bench/
  README.md
  run_three_way.py                    # full-pipeline tetra3 + LOST + ours
  run_decoupled_lost_vs_tetra3.py     # star-ID only: same centroids → both solvers
  export_centroids.py
  export_hipparcos_for_lost.py        # Hipparcos → LOST TSV catalog
  data/
    hipparcos_for_lost.tsv
    lost_database_hip.dat             # LOST k-vector DB built from Hipparcos
  centroids/
    lost_native/                      # LOST-generated frozen centroids + truth
  adapters/
    tetra3_from_centroids.py
    lost_from_centroids.cpp           # pyramid+DQM on frozen centroids
    lost_from_centroids               # WSL binary
    Makefile.lost_adapter
    lost_run.sh
  results/
    DECOUPLED_LOST_VS_TETRA3.md
    decoupled_lost_vs_tetra3.json
```

## Decoupled LOST vs tetra3 (fair star-ID comparison)

```powershell
# Build adapter (WSL)
wsl -e bash -lc "cd /mnt/c/Users/Reddy/Desktop/st/star_tracker_suite/bench/adapters && make -f Makefile.lost_adapter"

# Run on LOST-native frozen centroids (recommended)
python bench/run_decoupled_lost_vs_tetra3.py `
  --centroid-dir bench/centroids/lost_native `
  --extractor lost_generate --fov 20
```

Latest result (5 cases, identical centroids): both solvers **5/5**, **0 false solves**,
median attitude disagreement **~0.0004 deg**.

## Metrics (per frame)

- **solve** / **no-solve** / **false-solve** (attitude error > threshold but reported valid)
- attitude error (arcsec): cross-boresight and about-boresight; roll reported separately
- wall time: p50, p95, p99.9, max
- peak RSS and static database size
- for ours (future): which frame in a 3-frame window produced the lock

False-solve rate is the metric that matters most. Our current pair-voting matcher cannot
measure it because it has no notion of "solve."

## Run (Windows host + WSL for LOST)

```powershell
# 1. Build our C++ benchmark runner (Visual Studio)
cmake -S . -B build-vs -G "Visual Studio 17 2022" -A x64
cmake --build build-vs --config Release --target star_tracker_benchmark_runner

# 2. Build LOST inside WSL (see benchmarking/lost_baseline/BASELINE.md)

# 3. Run three-way comparison
python bench/run_three_way.py --synthetic-cases 20
```

## Honest claim scope

Do not pitch "faster/smaller than tetra3" — LOST already beats tetra3 on both in print.
The contribution target is **multi-frame temporal consistency**: run pyramid at relaxed
per-frame thresholds and recover confidence across a short frame history, with bounded-time
anytime behavior for hard real-time scheduling.
