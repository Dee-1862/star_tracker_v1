"""Do refusals arrive alone, or in bursts?

Every other measurement in this project treats each frame as an independent
draw, which answers "how often does it refuse" but not "how long is it blind".
Those are different questions: a 1-in-300 refusal rate is harmless if the
refusals are scattered and serious if all 300 arrive consecutively, and the two
are indistinguishable in any per-frame statistic.

This simulates a continuously slewing spacecraft, so successive frames are
correlated the way real ones are, and reports the distribution of CONSECUTIVE
refusals. Gyro coasting time, not refusal rate, is what an ADCS engineer needs.
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

EXE = str(ROOT / "build-vs" / "flight_software" / "Release" /
          "star_tracker_sequence_runner.exe")
W = H = 1024
FOV = 20.0
NGP_RA, NGP_DEC = 192.85948, 27.12825


def galactic_latitude(ra_deg: float, dec_deg: float) -> float:
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    rp, dp = math.radians(NGP_RA), math.radians(NGP_DEC)
    s = (math.sin(dec) * math.sin(dp) +
         math.cos(dec) * math.cos(dp) * math.cos(ra - rp))
    return math.degrees(math.asin(max(-1.0, min(1.0, s))))


def runs_of_zeros(flags: list[int]) -> list[int]:
    """Lengths of consecutive refusals."""
    out, cur = [], 0
    for f in flags:
        if f:
            if cur:
                out.append(cur)
            cur = 0
        else:
            cur += 1
    if cur:
        out.append(cur)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frames", type=int, default=400)
    ap.add_argument("--rate-hz", type=float, default=10.0)
    ap.add_argument("--slew-deg-s", type=float, default=0.30,
                    help="Continuous slew. Successive frames overlap, as in "
                         "flight, instead of being independent samples.")
    args = ap.parse_args()

    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))

    # A great circle sweep crosses the galactic plane twice, so the run covers
    # both crowded and empty sky rather than an average of the two.
    step = args.slew_deg_s / args.rate_hz
    tmpdir = Path(tempfile.mkdtemp(prefix="burst_"))
    paths, truth = [], []
    print(f"rendering {args.frames} frames, {step:.4f} deg apart "
          f"({args.slew_deg_s} deg/s at {args.rate_hz} Hz)...")

    for i in range(args.frames):
        ang = math.radians(i * step)
        # Great circle inclined to the equator, so declination varies too.
        ra = math.degrees(math.atan2(math.sin(ang) * math.cos(math.radians(60.0)),
                                     math.cos(ang))) % 360.0
        dec = math.degrees(math.asin(math.sin(ang) * math.sin(math.radians(60.0))))
        roll = (i * step * 0.5) % 360.0

        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) >= 4:
            img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                               noise_sigma=0.5, psf_sigma=1.0, seed=i)
            c = np.asarray(tetra3.get_centroids_from_image(
                Image.fromarray(img), sigma=3, max_returned=30))
        else:
            c = np.empty((0, 2))
        f = tmpdir / f"f{i:05d}.tsv"
        with f.open("w") as fh:
            fh.write("# x y intensity\n")
            for j, row in enumerate(c):
                fh.write(f"{row[1]:.6f} {row[0]:.6f} {3000 - j}\n")
        paths.append(f)
        truth.append((ra, dec, len(c), abs(galactic_latitude(ra, dec))))
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{args.frames}", flush=True)

    listfile = tmpdir / "list.txt"
    listfile.write_text("\n".join(str(p) for p in paths) + "\n", encoding="utf-8")

    out = subprocess.run([EXE, str(listfile), str(W), str(H), f"{FOV}"],
                         capture_output=True, text=True).stdout
    rows = [l.split(",") for l in out.strip().splitlines()[1:]]
    if not rows:
        print("no output from sequence runner")
        return

    solved = [int(r[1]) for r in rows]
    modes = [r[2] for r in rows]
    reasons = [r[3] for r in rows]
    tracked = [float(r[11]) for r in rows]

    n = len(solved)
    nsolved = sum(solved)
    bursts = runs_of_zeros(solved)
    period = 1.0 / args.rate_hz

    print(f"\n=== CONTINUOUS SWEEP: {n} frames at {args.rate_hz} Hz "
          f"({n * period:.0f} s of flight) ===\n")
    print(f"  solved            {nsolved}/{n}  ({100.0 * nsolved / n:.1f}%)")
    print(f"  refused           {n - nsolved}")
    print(f"  acquisition mode  {sum(1 for m in modes if m == 'acquisition')} frames")

    if bursts:
        b = np.array(bursts)
        print(f"\n  refusal bursts    {len(b)}")
        print(f"    longest         {b.max()} frames = {b.max() * period:.2f} s BLIND")
        print(f"    median          {int(np.median(b))} frames")
        print(f"    p95             {int(np.percentile(b, 95))} frames")
        hist = {}
        for x in b:
            hist[int(x)] = hist.get(int(x), 0) + 1
        print("    length: count   " +
              "  ".join(f"{k}:{v}" for k, v in sorted(hist.items())[:12]))
    else:
        print("\n  refusal bursts    none")

    if n - nsolved:
        rc = {}
        for s, r in zip(solved, reasons):
            if not s:
                rc[r] = rc.get(r, 0) + 1
        print(f"\n  refusal reasons   {rc}")
        stars = [t[2] for s, t in zip(solved, truth) if not s]
        gb = [t[3] for s, t in zip(solved, truth) if not s]
        if stars:
            print(f"  refused frames    median {int(np.median(stars))} centroids, "
                  f"median |b| {np.median(gb):.0f} deg")
        sok = [t[2] for s, t in zip(solved, truth) if s]
        if sok:
            print(f"  solved frames     median {int(np.median(sok))} centroids")

    tr = np.array([x for x in tracked if x > 0])
    if len(tr):
        print(f"\n  focal tracked     {tr[0]:.2f} -> {tr[-1]:.2f} px "
              f"(spread {tr.max() - tr.min():.2f} px)")

    print("\n  Independent-frame statistics cannot see any of this: they report")
    print("  the refusal RATE, while an ADCS needs the longest BLIND INTERVAL.")


if __name__ == "__main__":
    main()
