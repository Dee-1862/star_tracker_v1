# Three-Way Star Tracker Benchmark Report

## Solvers

- **Ours:** allocation-free C++ pair-voting matcher (Release x64)
- **tetra3:** ESA 4-star pattern hash + binomial verify
- **LOST:** pyramid + k-vector + DQM attitude (UWCubeSat)

## Results

### Synthetic (20 cases, 20 deg FOV, Hipparcos truth)

| Solver | Field solves | False solves | Median latency |
|--------|--------------|--------------|----------------|
| Ours   | 0/20 | n/a (86.75% false IDs) | 2.789 ms |
| tetra3 | 19/20 | 0 | 90.404 ms |
| LOST   | 20/20 | 20 | 461.706249974668 ms |

### Real images (2 ESA FLIR frames)

- tetra3 solved: 2/2
- Ours IDs agreeing with tetra3: 0 (disagreed: 40)
- LOST: see lost_summary in JSON

## LOST detail

- LOST synthetic solves: 20/20; false solves: 20; crashes: 0.
- LOST real-image solves: 0/2; median agreement vs tetra3: None deg.
- LOST median wall time (synthetic): 461.706249974668 ms; database size: 1125.9 KiB.

Raw per-case data: `three_way_results.json`
