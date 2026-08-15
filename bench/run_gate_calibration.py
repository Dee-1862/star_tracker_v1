"""Re-derive the integrity threshold now that the residual metric is sound.

Sweeps uniform sky plus four contamination modes and reports the residual
distribution of correct versus incorrect solves.
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
TRUE_THR = 0.05
FIELDS = 40
# Deliberately wide so the gate never fires: we want the raw residual of every
# solve, correct or not, in order to choose the threshold from the data.
OPEN_GATE = "100000"


def run(xy: np.ndarray, fov: float) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        tsv = Path(f.name)
        f.write("# x y intensity\n")
        for j, (x, y) in enumerate(xy):
            f.write(f"{x:.6f} {y:.6f} {2000 - j}\n")
    out = subprocess.run([EXE, str(tsv), str(W), str(H), f"{fov}", OPEN_GATE],
                         capture_output=True, text=True).stdout
    tsv.unlink(missing_ok=True)
    return {p[0]: p[1] for p in (l.split() for l in out.splitlines()) if len(p) == 2}


def main() -> None:
    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    rng = np.random.default_rng(2024)

    rows: list[tuple[str, str, float, int]] = []  # mode, verdict, residual, clique

    for i in range(FIELDS):
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 8:
            continue
        img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                           noise_sigma=0.5, psf_sigma=1.0, seed=9000 + i)
        c = np.asarray(tetra3.get_centroids_from_image(
            Image.fromarray(img), sigma=3, max_returned=30))
        if len(c) < 8:
            continue
        xy = np.column_stack((c[:, 1], c[:, 0]))

        cases: list[tuple[str, np.ndarray, float]] = [("clean", xy, FOV)]

        hot = np.vstack([xy, rng.uniform(20, W - 20, size=(6, 2))])
        cases.append(("hot_pixels", hot, FOV))

        for pct in (-5.0, -2.0, 2.0, 5.0):
            cases.append((f"fov{pct:+.0f}%", xy, FOV * (1.0 + pct / 100.0)))

        cases.append(("pure_noise", rng.uniform(20, W - 20, size=(25, 2)), FOV))

        for mode, points, fov in cases:
            d = run(points, fov)
            clique = int(d.get("clique_size", 0))
            if d.get("attitude_known") != "1":
                rows.append((mode, "REFUSE", -1.0, clique))
                continue
            err = angular_error_degrees(float(d["attitude_ra"]),
                                        float(d["attitude_de"]), ra, dec)
            verdict = "TRUE" if err < TRUE_THR else "FALSE"
            rows.append((mode, verdict, float(d["residual_rms_arcsec"]), clique))

    def stats(sel):
        v = np.array([r[2] for r in rows if sel(r) and r[2] >= 0.0])
        return v

    print(f"\nfields: {FIELDS}   samples: {len(rows)}\n")
    print(f"{'mode':>12} {'TRUE':>5} {'FALSE':>6} {'REFUSE':>7} "
          f"{'resid p50':>10} {'resid max':>10} {'clique p50':>11}")
    print("-" * 68)
    for mode in ("clean", "hot_pixels", "fov-5%", "fov-2%", "fov+2%", "fov+5%",
                 "pure_noise"):
        sub = [r for r in rows if r[0] == mode]
        if not sub:
            continue
        t = sum(1 for r in sub if r[1] == "TRUE")
        f_ = sum(1 for r in sub if r[1] == "FALSE")
        x = sum(1 for r in sub if r[1] == "REFUSE")
        res = np.array([r[2] for r in sub if r[2] >= 0.0])
        cl = np.array([r[3] for r in sub if r[1] != "REFUSE"])
        p50 = f"{np.median(res):.2f}" if len(res) else "-"
        mx = f"{res.max():.2f}" if len(res) else "-"
        c50 = f"{int(np.median(cl))}" if len(cl) else "-"
        print(f"{mode:>12} {t:>5} {f_:>6} {x:>7} {p50:>10} {mx:>10} {c50:>11}")

    true_r = stats(lambda r: r[1] == "TRUE")
    false_r = stats(lambda r: r[1] == "FALSE")

    print("\n--- residual distribution, all modes pooled (arcsec) ---")
    for name, v in (("TRUE", true_r), ("FALSE", false_r)):
        if len(v):
            print(f"  {name:>5}  n={len(v):<4} min {v.min():8.2f}  p50 {np.median(v):8.2f}"
                  f"  p95 {np.percentile(v, 95):8.2f}  max {v.max():8.2f}")
        else:
            print(f"  {name:>5}  n=0")

    if len(true_r) and len(false_r):
        hi, lo = true_r.max(), false_r.min()
        print(f"\n  worst TRUE  {hi:.2f}   best FALSE {lo:.2f}")
        if hi < lo:
            print(f"  SEPARABLE, margin {lo / hi:.1f}x -> threshold "
                  f"{math.sqrt(hi * lo):.2f} arcsec")
        else:
            print("  OVERLAP")
    elif len(true_r):
        print(f"\n  no false solves produced in any mode; worst TRUE residual "
              f"{true_r.max():.2f} arcsec")
        print(f"  suggested threshold (4x worst TRUE): {4 * true_r.max():.1f} arcsec")


if __name__ == "__main__":
    main()
