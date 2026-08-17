"""Maximum slew rate: how fast can the spacecraft rotate and still get a fix?

A rotating spacecraft smears each star into a streak during the exposure. Past
some rate the centroider rejects them as too elongated and the tracker refuses.
Every star tracker datasheet quotes this number; it is among the first things an
ADCS engineer asks for.

Blur is modelled by rendering each star as a series of sub-exposures along its
motion vector, which is what a real streak is.
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
FOCAL_PX = (W / 2.0) / math.tan(math.radians(FOV) / 2.0)
SUBSTEPS = 9


def smear(pixels, mags, rate_deg_s, exposure_s, angle_rad):
    """Split each star into SUBSTEPS dimmer copies along its motion vector."""
    travel_px = math.radians(rate_deg_s * exposure_s) * FOCAL_PX
    if travel_px < 0.05:
        return pixels, mags
    dx = math.cos(angle_rad) * travel_px
    dy = math.sin(angle_rad) * travel_px
    # Each sub-exposure carries 1/N of the light: +2.5*log10(N) magnitudes.
    dim = 2.5 * math.log10(SUBSTEPS)
    out_p, out_m = [], []
    for s in range(SUBSTEPS):
        f = (s / (SUBSTEPS - 1.0)) - 0.5
        out_p.append(pixels + np.array([dx * f, dy * f]))
        out_m.append(mags + dim)
    return np.vstack(out_p), np.concatenate(out_m)


def run(xy):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        tsv = Path(f.name)
        f.write("# x y intensity\n")
        for j, (x, y) in enumerate(xy):
            f.write(f"{x:.6f} {y:.6f} {3000 - j}\n")
    out = subprocess.run([EXE, str(tsv), str(W), str(H), f"{FOV}", "30", "0"],
                         capture_output=True, text=True).stdout
    tsv.unlink(missing_ok=True)
    return {p[0]: p[1] for p in (l.split() for l in out.splitlines()) if len(p) == 2}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields", type=int, default=20)
    ap.add_argument("--exposure", type=float, default=0.100)
    args = ap.parse_args()

    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    rng = np.random.default_rng(4242)

    scenes = []
    while len(scenes) < args.fields:
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 10:
            continue
        scenes.append((ra, dec, px, mags, float(rng.uniform(0, 2 * math.pi))))

    rates = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    print(f"\nexposure {args.exposure * 1000:.0f} ms, {len(scenes)} fields, "
          f"focal {FOCAL_PX:.0f} px")
    print(f"{'rate':>8} {'smear':>8} {'centroids':>10} {'solved':>7} "
          f"{'wrong':>6} {'refused':>8} {'err mdn':>9} {'clique':>7}")
    print("-" * 70)

    for rate in rates:
        travel = math.radians(rate * args.exposure) * FOCAL_PX
        solved = wrong = refused = 0
        errs, cents, cliques = [], [], []
        for ra, dec, px, mags, ang in scenes:
            p, m = smear(px, mags, rate, args.exposure, ang)
            img = render_image(p, m, width=W, height=H, noise_mean=2.0,
                               noise_sigma=0.5, psf_sigma=1.0, seed=7)
            c = np.asarray(tetra3.get_centroids_from_image(
                Image.fromarray(img), sigma=3, max_returned=30))
            cents.append(len(c))
            if len(c) < 5:
                refused += 1
                continue
            d = run(np.column_stack((c[:, 1], c[:, 0])))
            if d.get("attitude_known") != "1":
                refused += 1
                continue
            e = angular_error_degrees(float(d["attitude_ra"]),
                                      float(d["attitude_de"]), ra, dec)
            cliques.append(int(d.get("clique_size", 0)))
            if e < 0.05:
                solved += 1
                errs.append(e)
            else:
                wrong += 1
        em = f"{np.median(errs):.4f}" if errs else "-"
        cq = f"{int(np.median(cliques))}" if cliques else "-"
        print(f"{rate:>7.1f}°/s {travel:>7.1f}px {int(np.median(cents)):>10} "
              f"{solved:>7} {wrong:>6} {refused:>8} {em:>9} {cq:>7}")

    print(f"\nmax slew rate = highest rate with full solve rate and zero wrong "
          f"answers.\nScales inversely with exposure: halve the exposure, "
          f"double the rate.")


if __name__ == "__main__":
    main()
