# Star Tracker Benchmark Report

## Scope and provenance

- C++ candidate: this repository's allocation-free centroiding, three-frame kinematic filter, and pair-distance LIS matcher, Release x64 build.
- Baseline: ESA `tetra3` commit `f9fa2eb9a32a5efc529e2d86f0b59f35b1e9028d`.
- Synthetic validation: 20 deterministic 20-degree scenes rendered from the real ESA Hipparcos main catalogue (`hip_main.dat`, magnitude <= 6). These are simulated detector images, not real captures.
- Real-image validation: 2 official ESA tetra3 FLIR Blackfly test frames. This is a small smoke test, not a statistically representative flight qualification campaign.

## Results

- Synthetic field solves: C++ 0/20; tetra3 19/20.
- C++ star IDs: 53 correct, 347 incorrect; false-ID rate 86.75%.
- tetra3 false solves: 0; median pointing error 0.012856 deg.
- Median synthetic end-to-end host latency: C++ 2.789 ms (index construction excluded); tetra3 90.404 ms.
- Memory: C++ core objects plus catalogue 443.4 KiB; C++ host process peak 5.7 MiB. tetra3 database arrays 94.6 MiB; Python RSS after database load 266.0 MiB. The C++ host-process figure includes runtime/OS overhead and the 1 MiB input buffer; the core figure is the relevant static algorithm budget.
- Official real frames solved by tetra3: 2/2. C++ IDs agreeing with tetra3: 0; disagreeing: 40.

## Interpretation

The benchmark is now operational and reproducible, but the two systems are not equivalent. tetra3 returns and verifies a full attitude solution using a four-star pattern database. The current C++ pair-voting matcher returns per-centroid IDs without a geometric attitude verification stage. Any nonzero but incorrect C++ IDs are therefore false identifications, not successful lost-in-space solutions.

Real-frame "agreement" uses tetra3's matched Hipparcos IDs as the reference because the two bundled images do not include independent truth files. It demonstrates execution on real sensor data but must not be presented as independent accuracy validation.

Raw per-case measurements and configuration are in `comparison_results.json`.
