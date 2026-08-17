"""Statistical power run: bound the false-solve rate, don't just report zero.

Zero failures in n trials does not mean zero rate. The one-sided 95% upper
bound (Clopper-Pearson; the "rule of three" 3/n is its approximation) is the
number a datasheet can carry. This sweeps a wide operating envelope and reports
that bound, both overall and per condition.

Centroids are cached, so re-running with a larger --fields only renders the new
ones.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
for p in ("benchmarking", "sim_environment/src", "benchmarking/tetra3_baseline"):
    sys.path.insert(0, str(ROOT / p))
np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from generate_star_field import camera_basis, project_catalog, render_image  # noqa: E402
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402
from run_comparison import angular_error_degrees  # noqa: E402

EXE = str(ROOT / "build-vs" / "flight_software" / "Release" /
          "star_tracker_centroid_runner.exe")
CACHE = ROOT / "bench" / "data" / "false_solve_scenes.npz"
W = H = 1024
FOV = 20.0
POS_THR = 0.05


def clopper_pearson_upper(failures: int, trials: int, confidence: float = 0.95) -> float:
    """One-sided upper bound on the failure probability."""
    if trials == 0:
        return 1.0
    if failures == 0:
        return 1.0 - (1.0 - confidence) ** (1.0 / trials)
    return float(stats.beta.ppf(confidence, failures + 1, trials - failures))


def build_scenes(fields: int, rng) -> list:
    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    scenes = []
    attempts = 0
    while len(scenes) < fields and attempts < fields * 4:
        attempts += 1
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 8:
            continue
        img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                           noise_sigma=0.5, psf_sigma=1.0, seed=len(scenes) * 7 + 1)
        c = np.asarray(tetra3.get_centroids_from_image(
            Image.fromarray(img), sigma=3, max_returned=30))
        if len(c) < 8:
            continue
        scenes.append((ra, dec, np.column_stack((c[:, 1], c[:, 0]))))
        if len(scenes) % 25 == 0:
            print(f"  rendered {len(scenes)}/{fields}", flush=True)
    return scenes


def run(xy: np.ndarray, fov: float, search: float) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        tsv = Path(f.name)
        f.write("# x y intensity\n")
        for j, (x, y) in enumerate(xy):
            f.write(f"{x:.6f} {y:.6f} {3000 - j}\n")
    out = subprocess.run([EXE, str(tsv), str(W), str(H), f"{fov}", "30", f"{search}"],
                         capture_output=True, text=True).stdout
    tsv.unlink(missing_ok=True)
    return {p[0]: p[1] for p in (l.split() for l in out.splitlines()) if len(p) == 2}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields", type=int, default=300)
    ap.add_argument("--rebuild-cache", action="store_true")
    args = ap.parse_args()

    rng = np.random.default_rng(20260815)

    if CACHE.is_file() and not args.rebuild_cache:
        z = np.load(CACHE, allow_pickle=True)
        scenes = list(z["scenes"])
        print(f"loaded {len(scenes)} cached scenes from {CACHE.name}")
    else:
        scenes = []
    if len(scenes) < args.fields:
        print(f"rendering {args.fields - len(scenes)} more scenes...")
        scenes = scenes + build_scenes(args.fields - len(scenes), rng)
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(CACHE, scenes=np.array(scenes, dtype=object))
    scenes = scenes[:args.fields]

    def mirror(xy):
        p = xy.copy()
        p[:, 1] = (H - 1) - p[:, 1]
        return p

    def contaminate(xy, n):
        return np.vstack([xy, rng.uniform(20, W - 20, size=(n, 2))])

    conditions = [
        ("clean",              lambda xy: xy,                 FOV,        0.0),
        ("mirrored",           mirror,                        FOV,        0.0),
        ("focal -1%",          lambda xy: xy,                 FOV * 0.99, 0.0),
        ("focal +1%",          lambda xy: xy,                 FOV * 1.01, 0.0),
        ("focal -3%",          lambda xy: xy,                 FOV * 0.97, 0.0),
        ("focal +3%",          lambda xy: xy,                 FOV * 1.03, 0.0),
        ("focal -2% +search",  lambda xy: xy,                 FOV * 0.98, 6.0),
        ("focal +2% +search",  lambda xy: xy,                 FOV * 1.02, 6.0),
        ("mirrored +search",   mirror,                        FOV,        6.0),
        ("6 false stars",      lambda xy: contaminate(xy, 6), FOV,        0.0),
        ("12 false stars",     lambda xy: contaminate(xy, 12), FOV,       0.0),
        ("pure noise",         lambda xy: rng.uniform(20, W - 20, size=(25, 2)),
                                                              FOV,        0.0),
    ]

    print(f"\n{len(scenes)} fields x {len(conditions)} conditions "
          f"= {len(scenes) * len(conditions)} attempts\n")
    print(f"{'condition':>18} {'attempts':>9} {'correct':>8} {'WRONG':>6} "
          f"{'refused':>8} {'95% upper bound':>16}")
    print("-" * 72)

    tot_att = tot_wrong = 0
    for name, transform, fov, search in conditions:
        att = correct = wrong = refused = 0
        for ra, dec, xy in scenes:
            d = run(transform(xy), fov, search)
            att += 1
            if d.get("attitude_known") != "1":
                refused += 1
                continue
            err = angular_error_degrees(float(d["attitude_ra"]),
                                        float(d["attitude_de"]), ra, dec)
            if err < POS_THR:
                correct += 1
            else:
                wrong += 1
        tot_att += att
        tot_wrong += wrong
        ub = clopper_pearson_upper(wrong, att)
        print(f"{name:>18} {att:>9} {correct:>8} {wrong:>6} {refused:>8} "
              f"{ub:>15.2e}")

    print("-" * 72)
    ub = clopper_pearson_upper(tot_wrong, tot_att)
    print(f"{'POOLED':>18} {tot_att:>9} {'':>8} {tot_wrong:>6} {'':>8} "
          f"{ub:>15.2e}")
    print(f"\nfalse-solve rate <= {ub:.2e} per attempt, 95% confidence "
          f"({tot_wrong} failures in {tot_att} attempts)")
    for target in (1e-2, 1e-3, 1e-4):
        need = math.ceil(math.log(0.05) / math.log(1.0 - target))
        print(f"  to claim <= {target:.0e} needs {need:,} consecutive "
              f"failure-free attempts")


if __name__ == "__main__":
    main()
