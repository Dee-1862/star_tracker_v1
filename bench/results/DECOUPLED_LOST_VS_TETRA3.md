# Decoupled Star-ID Benchmark: tetra3 vs LOST

## Setup

- Same frozen centroids for both solvers (`extractor=lost_generate`)
- FOV estimate: 20.0 deg
- Cases: 5
- LOST path: pyramid + DQM via `lost_from_centroids`
- tetra3 path: `solve_from_centroids`

## Results

| Solver | Solves | False solves | Median pointing error |
|--------|--------|--------------|-----------------------|
| tetra3 | 5/5 | 0 | 0.0004146501821704124 deg |
| LOST   | 5/5 | 0 | 3.1945284701301985e-06 deg |

- Both solved: 5/5
- Median attitude disagreement (when both solved): 0.00041847190749640963 deg
- LOST database size: 2696.4 KiB

Raw data: `decoupled_lost_vs_tetra3.json`
