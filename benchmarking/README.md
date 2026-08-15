# Reproducible comparison

This directory compares the allocation-free C++ flight pipeline with:

- ESA [tetra3](https://github.com/esa/tetra3) (4-star pattern hash + binomial verify)
- [LOST](https://github.com/UWCubeSat/lost) (pyramid + k-vector + DQM attitude)

The three-way harness lives in `../bench/`. See `../bench/README.md` for the decoupled
centroid benchmark design and LOST build instructions (`lost_baseline/BASELINE.md`).

## Inputs

- `../sim_environment/data/raw/hip_main.dat`: real Hipparcos catalogue source.
- `tetra3_baseline/tetra3/data/default_database.npz`: ESA's versioned tetra3
  database built from Hipparcos.
- `tetra3_baseline/examples/data/*.tiff`: ESA's two real FLIR Blackfly test
  images.

The random-attitude benchmark images are synthetic detector renders generated
from the real Hipparcos catalogue. They must not be described as real imagery.

## Run

From the repository root in a Visual Studio developer environment:

```powershell
cmake -S . -B build-vs -G "Visual Studio 17 2022" -A x64
cmake --build build-vs --config Release --target star_tracker_benchmark_runner
python benchmarking/run_comparison.py --synthetic-cases 20
```

Results are written to `results/comparison_results.json` and
`results/REPORT.md`.

tetra3 revision `f9fa2eb9a32a5efc529e2d86f0b59f35b1e9028d` uses
`np.math`, which NumPy 2 removed. The benchmark applies a local compatibility
alias; the ESA source remains unmodified.

## Metric boundaries

- A C++ field is successful only when it has at least three correct IDs and no
  incorrect IDs. A nonzero ID is not automatically counted as a match.
- A tetra3 synthetic solve is false when its boresight differs from scene truth
  by more than 0.1 degree.
- Real-image C++ IDs are compared with tetra3's matched Hipparcos IDs. That is
  baseline agreement, not independent ground truth.
- C++ pipeline latency excludes one-time index construction. Raw index timing is
  retained per case.
- The C++ static footprint and host process peak are reported separately. The
  latter includes the process runtime and 1 MiB camera input buffer.
