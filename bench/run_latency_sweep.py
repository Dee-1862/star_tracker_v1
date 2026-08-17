"""Latency and clique worst-case characterisation.

The consistency-graph search is the headline integrity feature and the biggest
timing unknown: maximum clique is NP-hard in general, and dense galactic-plane
fields are the adversarial input. This measures the actual work distribution,
stratified by galactic latitude, and reports the tail rather than the average.

A hard-real-time claim needs the max, not the median.
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
W = H = 1024
FOV = 20.0

# North galactic pole, J2000.
NGP_RA, NGP_DEC = 192.85948, 27.12825


def galactic_latitude(ra_deg: float, dec_deg: float) -> float:
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    ra_p, dec_p = math.radians(NGP_RA), math.radians(NGP_DEC)
    s = (math.sin(dec) * math.sin(dec_p) +
         math.cos(dec) * math.cos(dec_p) * math.cos(ra - ra_p))
    return math.degrees(math.asin(max(-1.0, min(1.0, s))))


def run(xy: np.ndarray, search: float) -> tuple[dict, float]:
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        tsv = Path(f.name)
        f.write("# x y intensity\n")
        for j, (x, y) in enumerate(xy):
            f.write(f"{x:.6f} {y:.6f} {3000 - j}\n")
    t0 = time.perf_counter()
    out = subprocess.run([EXE, str(tsv), str(W), str(H), f"{FOV}", "30", f"{search}"],
                         capture_output=True, text=True).stdout
    wall = (time.perf_counter() - t0) * 1000.0
    tsv.unlink(missing_ok=True)
    d = {p[0]: p[1] for p in (l.split() for l in out.splitlines()) if len(p) == 2}
    return d, wall


def pct(v, q):
    return float(np.percentile(v, q)) if len(v) else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields", type=int, default=250)
    args = ap.parse_args()

    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    rng = np.random.default_rng(884422)

    rows = []
    made = 0
    while made < args.fields:
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 8:
            continue
        img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                           noise_sigma=0.5, psf_sigma=1.0, seed=made)
        c = np.asarray(tetra3.get_centroids_from_image(
            Image.fromarray(img), sigma=3, max_returned=30))
        if len(c) < 8:
            continue
        made += 1
        xy = np.column_stack((c[:, 1], c[:, 0]))
        gb = abs(galactic_latitude(ra, dec))

        for mode, search in (("tracking", 0.0), ("acquisition", 6.0)):
            d, wall = run(xy, search)
            rows.append({
                "mode": mode, "gb": gb, "wall": wall,
                "stars": len(px), "centroids": len(c),
                "nodes": int(d.get("node_count", 0)),
                "exp": int(d.get("expansions", 0)),
                "clique": int(d.get("clique_size", 0)),
                "trials": int(d.get("focal_trials", 1)),
                "solved": d.get("attitude_known") == "1",
            })
        if made % 50 == 0:
            print(f"  {made}/{args.fields}", flush=True)

    for mode in ("tracking", "acquisition"):
        sub = [r for r in rows if r["mode"] == mode]
        e = np.array([r["exp"] for r in sub], dtype=float)
        w = np.array([r["wall"] for r in sub], dtype=float)
        n = np.array([r["nodes"] for r in sub], dtype=float)
        t = np.array([r["trials"] for r in sub], dtype=float)
        print(f"\n=== {mode.upper()} MODE  (n={len(sub)}) ===")
        print(f"{'metric':>22} {'p50':>10} {'p95':>10} {'p99':>10} {'max':>10}"
              f" {'max/p50':>9}")
        for label, v in (("clique expansions", e), ("graph nodes", n),
                         ("focal trials", t), ("wall time ms*", w)):
            ratio = (v.max() / pct(v, 50)) if pct(v, 50) > 0 else float("nan")
            print(f"{label:>22} {pct(v,50):>10.1f} {pct(v,95):>10.1f} "
                  f"{pct(v,99):>10.1f} {v.max():>10.1f} {ratio:>9.1f}x")

    print("\n* wall time includes ~35 ms one-time index build and process "
          "start per invocation;\n  it is an upper bound on in-flight cost, "
          "not a flight figure.")

    print("\n=== EXPANSIONS BY GALACTIC LATITUDE (tracking mode) ===")
    print("  dense Milky Way fields are the adversarial input for clique search")
    print(f"{'|b| band':>14} {'n':>5} {'stars p50':>10} {'exp p50':>10} "
          f"{'exp max':>10} {'clique p50':>11}")
    bands = [(0, 15), (15, 30), (30, 60), (60, 90)]
    track = [r for r in rows if r["mode"] == "tracking"]
    for lo, hi in bands:
        sub = [r for r in track if lo <= r["gb"] < hi]
        if not sub:
            continue
        e = np.array([r["exp"] for r in sub], dtype=float)
        s = np.array([r["stars"] for r in sub], dtype=float)
        cl = np.array([r["clique"] for r in sub], dtype=float)
        print(f"{f'{lo}-{hi} deg':>14} {len(sub):>5} {pct(s,50):>10.0f} "
              f"{pct(e,50):>10.0f} {e.max():>10.0f} {pct(cl,50):>11.0f}")

    solved = sum(1 for r in track if r["solved"])
    print(f"\ntracking-mode solves: {solved}/{len(track)}")
    theoretical = 80 * 80 * 20
    worst = max(r["exp"] for r in rows)
    print(f"structural cap on expansions: {theoretical:,} "
          f"(kMaxNodes^2 x kMaxCliqueObserved)")
    print(f"worst observed:               {worst:,} "
          f"({100.0 * worst / theoretical:.1f}% of cap)")


if __name__ == "__main__":
    main()
