# Star Tracker Suite

An allocation-free star tracker for small satellites, benchmarked head-to-head against
the two best-known open-source solvers on identical inputs.

The distinguishing property is not speed or accuracy. It is **knowing when to refuse**.

---

## 1. What a star tracker does, in plain language

A satellite needs to know which way it is pointing. The most accurate way to find out is
to photograph the stars, because the stars are fixed and their positions are known to
extreme precision.

So the job is:

```
photo of some bright dots  ->  which dot is which star?  ->  therefore, which way am I facing?
```

The hard part is the middle step. A photo of a 20-degree patch of sky shows perhaps 30
dots. The catalogue contains thousands of stars. Nothing about a single dot identifies it
— every star looks like a white blob. The only usable information is **the angles between
dots**, because rotating the camera does not change the angle between two stars.

Every star tracker exploits this. They differ in how.

### The one term you need: focal length

Focal length is how zoomed-in the camera is. It converts *"these two dots are 200 pixels
apart"* into *"these two stars are 4 degrees apart in the sky"*. Get it wrong and every
angle you compute is wrong by the same percentage — a ruler with mis-printed markings.
This single number drives most of the differences below.

---

## 2. The three solvers

### tetra3 (ESA) — match the *shape*

Takes 4 dots, measures the 6 gaps between them, and **divides them all by the largest**:

```
gaps:    2.1°  3.4°  4.0°  5.2°  6.8°  8.5°
÷ 8.5 →  0.25  0.40  0.47  0.61  0.80      <- looked up in a dictionary of every 4-star shape
```

It stores no absolute angle, only proportions. If your focal length is wrong by 5%, every
gap shrinks by 5% and **the proportions are unchanged** — the error cancels. It can even
work the focal length out afterwards.

> Like recognising a face in a photo. You do not need to know whether the print is
> postcard or poster sized; the proportions are the same.

The price is the dictionary: **49.4 MB**, 6.18 million patterns. That does not fit on a
satellite.

### LOST (UWCubeSat) — measure the *distances*, then get a witness

Takes three dots, measures the **actual** angles, and finds a matching triangle in a sorted
list of every catalogue pair. Then brings in a fourth dot as a witness: if its three
distances to the first three also check out, the answer is accepted. If two different star
groups both fit, it discards the whole thing rather than guessing.

> Like identifying a city from GPS coordinates. Extremely precise — if your GPS is
> accurate. If it is off by even a little you confidently land in the wrong city, because
> wrong coordinates are still perfectly self-consistent coordinates.

Database: **2.7 MB**. But the ruler must be right.

### Ours — find the *largest group whose stories agree*

Instead of a fixed group of four, build a graph: each node is a guess ("dot 7 might be
star HIP 32349"), and two nodes are connected when they agree with each other. Then find
the **largest set in which every guess agrees with every other guess**.

Three consequences:

- **It scales to the sky.** A dense field yields a 17-star mutually-consistent group,
  which is overwhelming evidence. A sparse field yields 5, accepted with less confidence.
  Neither is a special case.
- **Outliers fall out for free.** A hot pixel agrees with nothing, so it never joins the
  group. Measured below at 50% contamination.
- **The group size is the confidence.** No separate score to calibrate.

Then two checks the others lack in combination:

1. **Handedness.** Flipping an image leaves every distance unchanged — hold a photo to a
   mirror and every gap between features is identical. Distances alone are blind to it.
   We check the sign of a triple product, which is not.
2. **Refusal.** After solving, predict where every matched star should appear and measure
   the disagreement. Above a calibrated threshold, report *no answer* rather than a wrong
   one.

### The pipeline

```
 [0] frame                 1024x1024, 8-bit
      |
 [1] centroiding           connected components, one raster pass, no recursion
      |                    -> ~30 dot positions
 [2] conditioning          brightest 20, 3-frame kinematic outlier filter
      |
 [3] bearing projection    pixel -> unit vector.  <-- focal length enters HERE, once,
      |                    explicitly. Both baselines bury this inside their matcher,
      |                    which is why a focal error there becomes a wrong answer
      |                    instead of a refusal.
      |
 [4] shortlist             pair-angle lookup -> up to 4 candidate stars per dot
      |
 [5] consistency graph     edge = two guesses agree; + handedness check
      |
 [6] largest clique        the mutually-agreeing set. size = confidence
      |
 [7] attitude              QUEST (no SVD, no eigensolver, no matrix library)
      |
 [8] integrity gate        reprojection residual > 30 arcsec -> REFUSE
      |
 [9] output                quaternion, RA/Dec/Roll, residual, refusal reason
```

Stage 3 is the structural difference from both baselines.

---

## 3. Results

All numbers below are reproducible with the commands in section 5. Both baselines are
**unmodified** — see section 4.

### 3.1 Head to head

240 sky positions (HEALPix nside=2) x 5 roll angles = 480 cases per solver.
Rendered star fields, identical extracted centroids fed to all three.
Correct = within 0.05 degrees.

| Input | Solver | Solved | Boresight correct | **Roll correct** |
|---|---|---|---|---|
| clean | **ours** | 240 | 240 | **240** |
| clean | tetra3 | 235 | 235 | 235 |
| clean | LOST | 240 | 240 | 240 |
| mirrored | **ours** | **0** | — | — |
| mirrored | tetra3 | 235 | 235 | **47** |
| mirrored | LOST | 240 | **0** | 1 |

On clean input all three are equivalent. On mirrored input:

- **LOST** reports success on all 240 and is wrong about where it is pointing on all 240
  (median error 92.9 degrees). It has a handedness check, and the check works — it
  correctly rejects the true match. But a failed check makes Pyramid *keep searching*, not
  stop, so it hunts until something else fits.
- **tetra3** gets the boresight right but the **roll wrong on 188 of 235**. It does not
  detect the flip; it absorbs it into roll. Its descriptor is reflection-invariant and its
  rotation fit has no determinant correction, so it returns a mirrored rotation that is
  perfectly self-consistent.
- **Ours refuses all 240.**

Neither baseline detects a mirrored image. LOST corrupts the boresight loudly; tetra3
corrupts the roll quietly. **This is only visible if you score roll** — scoring the
boresight alone makes tetra3 look flawless.

### 3.2 Accuracy (clean input)

| | Boresight error (median) | Roll error (median) |
|---|---|---|
| **ours** | 0.0128° | **11.3 arcsec** |
| tetra3 | 0.0130° | 19.6 arcsec |
| LOST | 0.0131° | (negated convention) |

Typical mutually-consistent group size: **13–17 stars**, against a required minimum of 5.

### 3.3 Robustness (ours)

40 sky positions, 280 samples:

| Condition | Correct | **Wrong** | Refused |
|---|---|---|---|
| clean | 40 | **0** | 0 |
| focal length −5% | 0 | **0** | 40 |
| focal length −2% | 0 | **0** | 40 |
| focal length +2% | 0 | **0** | 40 |
| focal length +5% | 0 | **0** | 40 |
| pure noise (25 random points) | 0 | **0** | 40 |

**Zero wrong answers in every condition.** Under focal error it refuses rather than
guessing. Fed random points it never invents a solve.

Focal sensitivity is a real limitation shared with LOST — absolute-angle matching needs a
correct ruler. tetra3 tolerates it because it uses proportions. See section 6.

### 3.4 Outlier tolerance (ours)

30 fields; injected false stars are made **brighter than every real star**, so they cannot
be filtered by brightness:

| False stars injected | Correct | **Wrong** | Refused |
|---|---|---|---|
| 0 | 30 | **0** | 0 |
| 3 | 30 | **0** | 0 |
| 6 | 30 | **0** | 0 |
| 10 (≈50% of stars considered) | 26 | **0** | 4 |

At 50% contamination it still solves 26/30 correctly and degrades into honest refusals.

### 3.5 Integrity threshold

Derived from data, not tuned. Over **190 correct solves** spanning uniform sky and up to
50% contamination, reprojection residual ranged **3.0 to 14.7 arcsec** (median 4.8).
Threshold set to **30 arcsec** — 2x the worst observed case.

The residual gate is a *backstop*, not the primary filter: the consistency and handedness
checks remove wrong answers upstream, and **no wrong answer reached the gate in any
sweep**.

### 3.6 Footprint

Measured with `sizeof`, x86-64 release:

| Component | Bytes |
|---|---|
| Matcher (pair index, graph, vote table) | 276,840 |
| Hipparcos catalogue, 5041 stars | 100,820 |
| Centroiding working state | 106,512 |
| Attitude solver | 24 |
| **Total static working state** | **484,196 (472.8 KB)** |

Plus a 1,048,576-byte frame buffer if the whole image is held in RAM. **The 500 KB budget
in `flight_software/LLM_SYSTEM_PROMPT.md` refers to working state and does not include the
frame** — a full-frame buffer alone exceeds it, so a flight build must stream or tile the
sensor read.

| | ours | tetra3 | LOST |
|---|---|---|---|
| Database | **101 KB** | 49.4 MB | 2.7 MB |
| Heap allocation | **none** | Python/NumPy | `std::vector` throughout |

---

## 4. Language, and what is actually ours

**This is C++11, not C.** If you need C specifically, say so — it is a mechanical but
non-trivial port (classes to structs, `std::array` to raw arrays, namespaces to prefixes).

Flight code obeys `flight_software/LLM_SYSTEM_PROMPT.md`. Audited feature use across all
of `flight_software/src` and `flight_software/include`:

| Feature | Count |
|---|---|
| `std::array` | 79 |
| `namespace` | 23 |
| `class` | 4 |
| `std::sort` (index sort at init, non-allocating) | 1 |
| `new` / `malloc` / `std::vector` / `std::string` | **0** |

Host-side tools (`benchmark_runner.cpp`, `centroid_runner.cpp`) do use `iostream` and
`std::string`. They never run on the satellite.

**Written here:**

| File | Purpose |
|---|---|
| `flight_software/src/centroiding.cpp` | connected-component centroiding |
| `flight_software/src/lis_grid_matcher.cpp` | pair index, consistency graph, clique search |
| `flight_software/src/attitude_solver.cpp` | bearing projection, QUEST, integrity gate |
| `flight_software/src/kinematic_outlier_filter.cpp` | 3-frame outlier rejection |
| `flight_software/centroid_runner.cpp` | frozen-centroid adapter (host tool) |
| `bench/run_*.py` | benchmark harnesses |

**Not ours** — cloned, unmodified, git-ignored: `benchmarking/lost_baseline/`,
`benchmarking/tetra3_baseline/`.

---

## 5. Reproducing the numbers

### 5.1 Prerequisites

- CMake ≥ 3.16 and a C++11 compiler (tested: MSVC 19.4x / Visual Studio 2022)
- Python 3.10+ with `numpy`, `pillow`, `scipy`
- WSL with `g++-12`, `libcairo2-dev`, `libeigen3-dev` (only to build LOST)

### 5.2 Fetch the baselines

Not committed. Pinned to the commits these results were produced against:

```bash
cd benchmarking

git clone https://github.com/UWCubeSat/lost lost_baseline
git -C lost_baseline checkout 6543ec8

git clone https://github.com/esa/tetra3 tetra3_baseline
git -C tetra3_baseline checkout f9fa2eb
```

tetra3 downloads its 49.4 MB database on first use. Note that `tetra3.py` calls
`np.math.factorial`, removed in NumPy ≥ 2.0; our adapters shim it with `np.math = math`.

### 5.3 External data

None of this is committed — together it is roughly **1.3 GB**, and git can never forget a
large file once added. Fetch what you need; the table says which results depend on each.

| Path (git-ignored) | Size | Needed for |
|---|---|---|
| `sim_environment/data/raw/hip_main.dat` | 51 MB | **Everything.** Hipparcos catalogue. |
| `sim_environment/data/raw/star-tracker-data.zip` | 20 MB | Real-image tests only |
| `benchmarking/dust/DUST.zip` | **822 MB** | Not used by any result below |
| `benchmarking/dust/DUST-code.zip` | 2.4 MB | Not used by any result below |
| `benchmarking/tetra3_baseline/tetra3/data/default_database.npz` | 49.4 MB | tetra3 comparison (auto-downloads) |
| `bench/data/lost_database_hip.dat` | 2.7 MB | LOST comparison (generated, §5.4) |

**Hipparcos** — the only mandatory download. The `hip_main.dat` catalogue from the ESA
Hipparcos archive (VizieR catalogue `I/239/hip_main`). Place it at
`sim_environment/data/raw/hip_main.dat`, then run §5.4 to build the processed binary, the
C++ static array, and the LOST database. Every number in section 3 comes from this.

**DUST** — *"DUST: An On-Orbit Star Tracker Benchmark for RSO Detection and Attitude
Estimation"*. Real on-orbit star tracker frames with pixel- and object-level annotations of
resident space objects (satellites and debris crossing the field), plus background
statistics and an annotation toolchain. `DUST.zip` is the imagery; `DUST-code.zip` is the
processing and annotation code.

> **Nothing in this README uses DUST yet.** It is present in the working tree as a
> candidate for future work and no result depends on it. It is also, at 822 MB, by far
> the largest thing here.

It is the natural next validation step: section 3.4 measures outlier tolerance against
*synthetic* injected false stars, and DUST provides real ones with ground-truth
annotations. Section 6 lists "synthetic imagery only" as a limitation; DUST is how that
gets retired.

> **Source not recorded.** The archives contain no URL, DOI, or citation file, and nothing
> in this repository documents where they were obtained. If you are the person who
> downloaded them, please replace this note with the dataset URL and citation before
> publishing — results should not depend on a dataset whose provenance we cannot state.

Extract to `benchmarking/dust/` (the whole directory is git-ignored):

```bash
mkdir -p benchmarking/dust
# place DUST.zip and DUST-code.zip here, then:
cd benchmarking/dust && unzip -q DUST.zip && unzip -q DUST-code.zip
```

### 5.4 Build

```powershell
# Ours
cmake -S . -B build-vs -G "Visual Studio 17 2022" -A x64
cmake --build build-vs --config Release

# Tests (4/4 expected)
ctest --test-dir build-vs -C Release --output-on-failure
```

```bash
# LOST, inside WSL
wsl -e bash -lc "cd /mnt/c/<path>/star_tracker_suite/benchmarking/lost_baseline && make -j"
wsl -e bash -lc "cd /mnt/c/<path>/star_tracker_suite/bench/adapters && make -f Makefile.lost_adapter"
```

### 5.5 Catalogue and database

```bash
# Hipparcos -> processed binary -> C++ static array
python sim_environment/src/prepare_catalog.py
python sim_environment/src/generate_cpp_catalog.py

# LOST k-vector database from the same Hipparcos stars
python bench/export_hipparcos_for_lost.py
```

### 5.6 Run the benchmarks

```bash
# Section 3.1 / 3.2 -- head to head, scores BOTH boresight and roll
python bench/run_roll_sweep.py

# Quick 12-field version of the same comparison
python bench/run_three_way_decoupled.py

# Section 3.3 / 3.5 -- robustness sweep and threshold derivation
python bench/run_gate_calibration.py
```

Single frame through our solver:

```bash
build-vs/flight_software/Release/star_tracker_centroid_runner.exe \
    centroids.tsv 1024 1024 20.0 [MAX_RESIDUAL_ARCSEC=30]
```

`centroids.tsv` is `x y intensity` per line, top-left origin. Output is `key value` pairs
including `clique_size`, `residual_rms_arcsec`, `gate_reason`, `attitude_known`.

---

## 6. Honest limitations

- **Focal sensitivity.** Refuses under ±2% focal error, exactly like LOST. tetra3 tolerates
  ±5% because proportions cancel scale. The intended fix is a two-path design: a
  scale-invariant method for occasional calibration, ours for every frame afterwards. Not
  built.
- **The integrity gate is unfalsified, not validated.** No wrong answer was produced in any
  sweep, so there is no measured false-solve rate — and therefore no basis for an
  aerospace-style integrity risk figure.
- **Synthetic imagery only.** Rendered fields with Gaussian PSF and Gaussian noise. No real
  sensor data, no optical distortion, no stray light, no motion blur.
- **Coordinate conventions are load-bearing and were wrong twice.** The matcher's bearing
  frame was left-handed relative to the catalogue (determinant −1), invisible because
  dot products ignore handedness. And `generate_star_field.camera_basis` is 180° offset in
  roll from every solver — **any earlier result quoting roll error against this harness is
  suspect**.
- **`attitude_solver_test` was self-referential.** It built its truth with the same
  convention it tested, and passed throughout a bug that made the residual metric return
  ~65 arcsec of pure floating-point noise for any input. A regression case now checks the
  residual against a known injected perturbation, but the general risk stands.
- **Timing is not characterised.** No worst-case execution time, no percentile latency, no
  target-hardware measurement. Nothing here supports a real-time scheduling claim.

---

## 7. Layout

```
flight_software/       constrained C++11 -- no heap, fixed arrays
  include/star_tracker/
  src/
  test/                4 test executables, run via ctest
  centroid_runner.cpp  frozen-centroid adapter (host)
  benchmark_runner.cpp full-image runner (host)

bench/                 comparison harnesses
  run_roll_sweep.py            head to head, both axes
  run_three_way_decoupled.py   12-field quick comparison
  run_gate_calibration.py      robustness + threshold derivation
  adapters/                    LOST frozen-centroid adapter
  data/                        generated databases (git-ignored)

sim_environment/       catalogue prep and synthetic star field rendering
benchmarking/          cloned baselines (git-ignored)
```
