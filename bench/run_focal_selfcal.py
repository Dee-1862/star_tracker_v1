"""Focal-length self-calibration: can ours recover from a lens that changed?

Simulates thermal drift by telling the solver a focal length that is wrong by a
known percentage, then checks whether it (a) refuses, (b) reacquires by search,
and (c) recovers the true focal length from the matched stars.
"""
from __future__ import annotations

import math
import subprocess
import sys
import tempfile
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
from run_comparison import angular_error_degrees  # noqa: E402

EXE = str(ROOT / "build-vs" / "flight_software" / "Release" /
          "star_tracker_centroid_runner.exe")
W = H = 1024
FOV = 20.0
TRUE_FOCAL = (W / 2.0) / math.tan(math.radians(FOV) / 2.0)
DRIFTS = (0.0, -0.5, -1.0, -2.0, -3.0, 1.0, 2.0, 3.0, 5.0)
FIELDS = 12
SEARCH_PCT = 6.0


def run(xy: np.ndarray, fov: float, search: float) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        tsv = Path(f.name)
        f.write("# x y intensity\n")
        for j, (x, y) in enumerate(xy):
            f.write(f"{x:.6f} {y:.6f} {3000 - j}\n")
    out = subprocess.run(
        [EXE, str(tsv), str(W), str(H), f"{fov}", "30", f"{search}"],
        capture_output=True, text=True).stdout
    tsv.unlink(missing_ok=True)
    return {p[0]: p[1] for p in (l.split() for l in out.splitlines()) if len(p) == 2}


def main() -> None:
    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    rng = np.random.default_rng(31337)

    scenes = []
    while len(scenes) < FIELDS:
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 8:
            continue
        img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                           noise_sigma=0.5, psf_sigma=1.0, seed=6000 + len(scenes))
        c = np.asarray(tetra3.get_centroids_from_image(
            Image.fromarray(img), sigma=3, max_returned=30))
        if len(c) < 8:
            continue
        scenes.append((ra, dec, np.column_stack((c[:, 1], c[:, 0]))))

    print(f"\ntrue focal length: {TRUE_FOCAL:.2f} px   ({FIELDS} fields per row)\n")
    for search, label in ((0.0, "SEARCH OFF (current behaviour)"),
                          (SEARCH_PCT, f"SEARCH ON (+/-{SEARCH_PCT:.0f}%)")):
        print(f"--- {label} ---")
        print(f"{'lens drift':>11} {'solved':>7} {'wrong':>6} {'refused':>8} "
              f"{'recovered focal':>16} {'error':>8} {'trials':>7}")
        for drift in DRIFTS:
            # A drifted lens means the TRUE focal length changed; the solver is
            # still configured with the nominal FOV, so it is wrong by -drift.
            told_fov = FOV * (1.0 + drift / 100.0)
            solved = wrong = refused = 0
            focals, trials = [], []
            for ra, dec, xy in scenes:
                d = run(xy, told_fov, search)
                trials.append(int(d.get("focal_trials", 1)))
                if d.get("attitude_known") != "1":
                    refused += 1
                    continue
                err = angular_error_degrees(float(d["attitude_ra"]),
                                            float(d["attitude_de"]), ra, dec)
                if err < 0.05:
                    solved += 1
                else:
                    wrong += 1
                focals.append(float(d["focal_refined_px"]))
            if focals:
                fm = float(np.median(focals))
                fe = 100.0 * (fm - TRUE_FOCAL) / TRUE_FOCAL
                fs, es = f"{fm:.2f}", f"{fe:+.3f}%"
            else:
                fs, es = "-", "-"
            print(f"{drift:>10.1f}% {solved:>7} {wrong:>6} {refused:>8} "
                  f"{fs:>16} {es:>8} {int(np.median(trials)):>7}")
        print()


if __name__ == "__main__":
    main()
