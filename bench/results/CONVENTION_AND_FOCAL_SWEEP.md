# Convention Diagnosis + Cross-Generator + Focal Sweep

## 1. Residual vectors on LOST-native centroids

RA/Dec residuals are small and change sign/magnitude across cases -> not a fixed convention bias; consistent with estimator/FOV-fit noise.

tetra3_roll ~ (-lost_roll) mod 360 on every case -> systematic roll-sign convention mismatch, not accuracy.

Per-case tetra3 residuals vs truth (deg):

| case | dRA | dDec | dRoll | dRoll vs -truth |
|------|-----|------|-------|-----------------|
| lost_synth_0000 | 0.000440 | 0.000030 | 0.0043 | 0.004347 |
| lost_synth_0001 | 0.000236 | 0.000212 | -80.0001 | -0.000051 |
| lost_synth_0002 | -0.000050 | -0.000009 | -160.0007 | -0.000749 |
| lost_synth_0003 | -0.004675 | -0.000286 | 119.9964 | -0.003583 |
| lost_synth_0004 | -0.000208 | -0.000499 | 40.0015 | 0.001485 |

- Mean tetra3 dRA/dDec: -0.000851 / -0.000110
- Mean tetra3−lost roll delta: -16.000°
- After negating LOST roll: 0.000290° (max abs 0.004362°)

## 2. Cross-generator (Hipparcos projector centroids)

### extractor=`perfect_project`
- tetra3: 5/5 (false 0), median err 0.0006950734802185297
- LOST: 5/5 (false 0), median err 0.0

### extractor=`tetra3_extract`
- tetra3: 5/5 (false 0), median err 0.012989938932367731
- LOST: 5/5 (false 0), median err 0.012989165308978775

## 3. Focal-length / FOV error sweep

Extractor: `perfect_project`, nominal FOV 20.0°

| FOV error | FOV used | tetra3 solves | tetra3 false | LOST solves | LOST false |
|-----------|----------|---------------|--------------|-------------|------------|
| -5.0% | 19.000 | 2/5 | 0 | 5/5 | 5 |
| -3.0% | 19.400 | 5/5 | 0 | 5/5 | 5 |
| -2.0% | 19.600 | 5/5 | 0 | 5/5 | 5 |
| -1.0% | 19.800 | 5/5 | 0 | 5/5 | 5 |
| +0.0% | 20.000 | 5/5 | 0 | 5/5 | 0 |
| +1.0% | 20.200 | 5/5 | 0 | 5/5 | 5 |
| +2.0% | 20.400 | 5/5 | 0 | 5/5 | 5 |
| +3.0% | 20.600 | 5/5 | 0 | 5/5 | 5 |
| +5.0% | 21.000 | 2/5 | 0 | 5/5 | 5 |
