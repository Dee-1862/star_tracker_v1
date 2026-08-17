# Calibration robustness sweep

- HEALPix nside=2, rolls=1, N~2
- Primary outcome threshold: 0.05° (sensitivity [0.01, 0.05, 0.2])
- tetra3 `fov_max_error` swept: [0.5, 1.0, 2.0]
- LOST oracle gates: False

tetra3 fov_max_error is swept; a cliff at ±(fov_max_error/FOV) is a search-window effect, not an intrinsic algorithm property. LOST rates use best-of-gates oracle when lost_gates=true.

## tetra3 TRUE/FALSE/REFUSE (primary thr), focal x fov_max_error

| focal % | fov_max=0.5° | fov_max=1° | fov_max=2° |
|---------|------|------|------|
| -5 | T0/F0/R2 | T1/F0/R1 | T2/F0/R0 |
| +0 | T2/F0/R0 | T2/F0/R0 | T2/F0/R0 |
| +5 | T0/F0/R2 | T0/F0/R2 | T2/F0/R0 |

## tetra3 timing (ms) at focal=0, vs fov_max_error

| fov_max_error | p50 | p95 | p99 | max |
|---------------|-----|-----|-----|-----|
| 0.5 | 3.5932500031776726 | 3.864914976293221 | 3.889062973903492 | 3.89509997330606 |
| 1 | 2.8033499838784337 | 2.805195050314069 | 2.8053590562194586 | 2.805400057695806 |
| 2 | 2.636949997395277 | 2.8898049844428897 | 2.9122809832915664 | 2.9178999830037355 |

## LOST best-of-gates (primary thr)

| focal % | TRUE | FALSE | REFUSE | true_rate |
|---------|------|-------|--------|-----------|
| -5 | 0 | 2 | 0 | 0.0 |
| +0 | 2 | 0 | 0 | 1.0 |
| +5 | 0 | 2 | 0 | 0.0 |
