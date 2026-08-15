# Calibration robustness sweep

- HEALPix nside=2, rolls=1, N≈12
- Primary outcome threshold: 0.05° (sensitivity [0.01, 0.05, 0.2])
- tetra3 `fov_max_error` swept: [0.5, 1.0, 2.0, 4.0]
- LOST oracle gates: True

tetra3 fov_max_error is swept; a cliff at ±(fov_max_error/FOV) is a search-window effect, not an intrinsic algorithm property. LOST rates use best-of-gates oracle when lost_gates=true.

## tetra3 TRUE/FALSE/REFUSE (primary thr) — focal × fov_max_error

| focal % | fov_max=0.5° | fov_max=1° | fov_max=2° | fov_max=4° |
|---------|------|------|------|------|
| -5 | T0/F0/R12 | T5/F0/R7 | T12/F0/R0 | T12/F0/R0 |
| -2 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| -1 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| -0.5 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| +0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| +0.5 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| +1 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| +2 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 | T12/F0/R0 |
| +5 | T0/F0/R12 | T8/F0/R4 | T12/F0/R0 | T12/F0/R0 |

## tetra3 timing (ms) at focal=0 — vs fov_max_error

| fov_max_error | p50 | p95 | p99 | max |
|---------------|-----|-----|-----|-----|
| 0.5 | 3.0934999813325703 | 10.821620066417381 | 10.857524083694443 | 10.866500088013709 |
| 1 | 3.740950021892786 | 10.293540061684325 | 10.79302802332677 | 10.91790001373738 |
| 2 | 3.2836999744176865 | 9.74476002738811 | 9.921111987205222 | 9.9651999771595 |
| 4 | 3.194499993696809 | 11.547074944246555 | 15.150454977992927 | 16.051299986429513 |

## LOST best-of-gates (primary thr)

| focal % | TRUE | FALSE | REFUSE | true_rate |
|---------|------|-------|--------|-----------|
| -5 | 0 | 0 | 12 | 0.0 |
| -2 | 0 | 0 | 12 | 0.0 |
| -1 | 8 | 0 | 4 | 0.6666666666666666 |
| -0.5 | 12 | 0 | 0 | 1.0 |
| +0 | 12 | 0 | 0 | 1.0 |
| +0.5 | 12 | 0 | 0 | 1.0 |
| +1 | 11 | 0 | 1 | 0.9166666666666666 |
| +2 | 0 | 0 | 12 | 0.0 |
| +5 | 0 | 0 | 12 | 0.0 |
