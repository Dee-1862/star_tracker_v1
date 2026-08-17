"""Speed and compute comparison: ours vs tetra3 vs LOST, identical centroids.

Uses each solver's OWN internal timer, so process startup, interpreter warmup
and one-time database loading are excluded. Those dominate naive wall-clock
measurement and say nothing about in-flight cost.

Startup cost is reported separately, because it is charged once and matters for
a different reason (boot-to-first-fix, not frame rate).
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
for p in ("benchmarking", "sim_environment/src", "benchmarking/tetra3_baseline"):
    sys.path.insert(0, str(ROOT / p))
np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from generate_star_field import camera_basis, project_catalog, render_image  # noqa: E402
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402

EXE = str(ROOT / "build-vs" / "flight_software" / "Release" /
          "star_tracker_centroid_runner.exe")
LOST = ROOT / "bench" / "adapters" / "lost_from_centroids"
LOSTDB = ROOT / "bench" / "data" / "lost_database_hip.dat"
W = H = 1024
FOV = 20.0


def wsl(p: Path) -> str:
    r = p.resolve()
    return f"/mnt/{r.drive.rstrip(':').lower()}{r.as_posix().split(':', 1)[1]}"


def kv(o: str) -> dict:
    return {p[0]: p[1] for p in (l.split() for l in o.splitlines()) if len(p) == 2}


def write(xy: np.ndarray) -> Path:
    f = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    f.write("# x y intensity\n")
    for j, (x, y) in enumerate(xy):
        f.write(f"{x:.6f} {y:.6f} {3000 - j}\n")
    f.close()
    return Path(f.name)


def pct(v, q):
    return float(np.percentile(v, q)) if len(v) else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields", type=int, default=60)
    args = ap.parse_args()

    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    rng = np.random.default_rng(99)

    scenes = []
    while len(scenes) < args.fields:
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 8:
            continue
        img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                           noise_sigma=0.5, psf_sigma=1.0, seed=len(scenes))
        c = np.asarray(tetra3.get_centroids_from_image(
            Image.fromarray(img), sigma=3, max_returned=30))
        if len(c) < 8:
            continue
        scenes.append(np.column_stack((c[:, 1], c[:, 0])))

    ours, lost, t3 = [], [], []
    ours_index, t3_load = [], []

    # tetra3 database load is one-time; measure it once, then reuse the object
    # so per-solve timings are not polluted by a 49.4 MB npz read.
    t0 = time.perf_counter()
    solver = tetra3.Tetra3()
    t3_load.append((time.perf_counter() - t0) * 1000.0)

    for xy in scenes:
        tsv = write(xy)

        d = kv(subprocess.run([EXE, str(tsv), str(W), str(H), f"{FOV}", "30", "0"],
                              capture_output=True, text=True).stdout)
        if "solve_ns" in d:
            ours.append(float(d["solve_ns"]) / 1e6)
        if "index_ns" in d:
            ours_index.append(float(d["index_ns"]) / 1e6)

        cmd = (f"{wsl(LOST)} {wsl(tsv)} {wsl(LOSTDB)} {W} {H} {FOV} 0.05 0.0001")
        dl = kv(subprocess.run(["wsl", "-e", "bash", "-lc", cmd],
                               capture_output=True, text=True).stdout)
        if "total_average_ns" in dl:
            lost.append(float(dl["total_average_ns"]) / 1e6)

        s = solver.solve_from_centroids(
            np.column_stack((xy[:, 1], xy[:, 0])), (H, W),
            fov_estimate=FOV, fov_max_error=1.0)
        if s.get("T_solve") is not None:
            t3.append(float(s["T_solve"]))

        tsv.unlink(missing_ok=True)

    print(f"\n{len(scenes)} identical centroid sets. Each solver's OWN timer.\n")
    print(f"{'solver':>10} {'n':>5} {'p50 ms':>9} {'p95 ms':>9} {'p99 ms':>9} "
          f"{'max ms':>9}")
    print("-" * 56)
    for name, v in (("ours", ours), ("tetra3", t3), ("LOST", lost)):
        a = np.array(v, dtype=float)
        if not len(a):
            print(f"{name:>10} {0:>5}   (unavailable)")
            continue
        print(f"{name:>10} {len(a):>5} {pct(a,50):>9.2f} {pct(a,95):>9.2f} "
              f"{pct(a,99):>9.2f} {a.max():>9.2f}")

    print("\n--- one-time startup (charged once, not per frame) ---")
    if ours_index:
        a = np.array(ours_index)
        print(f"{'ours':>10} catalogue pair index build: {a.mean():8.1f} ms")
    if t3_load:
        print(f"{'tetra3':>10} 49.4 MB database load:      {t3_load[0]:8.1f} ms")
    # Verified in lost_from_centroids.cpp: LoadDatabase and DeserializeCatalog
    # both complete BEFORE starIdBegin, so LOST's figure above is pure solve
    # and the comparison is like-for-like.
    print(f"{'LOST':>10} 2.7 MB database load:       (excluded, as for ours)")
    print("\n  CAVEAT: ours and tetra3 run natively; LOST runs under WSL2, which")
    print("  adds some overhead to its figures. Same x86-64 host throughout --")
    print("  none of this is target-hardware timing.")

    print("\n--- what the per-frame work actually is ---")
    print("  ours   : ~30k integer/float compatibility tests + one 3x3 QUEST")
    print("  tetra3 : ~32 hash probes x 70 patterns over a 12.4M-row table,")
    print("           then SVD + binomial test, in NumPy")
    print("  LOST   : 3 k-vector range queries + hash-map intersection,")
    print("           then Eigen 4x4 eigendecomposition (DQM)")


if __name__ == "__main__":
    main()
