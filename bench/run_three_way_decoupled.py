"""Decoupled three-way: ours vs tetra3 vs LOST on identical frozen centroids."""
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

OURS = ROOT / "build-vs" / "flight_software" / "Release" / "star_tracker_centroid_runner.exe"
LOST = ROOT / "bench" / "adapters" / "lost_from_centroids"
LOSTDB = ROOT / "bench" / "data" / "lost_database_hip.dat"
W = H = 1024
FOV = 20.0
THR = 0.05


def wsl(p: Path) -> str:
    r = p.resolve()
    return f"/mnt/{r.drive.rstrip(':').lower()}{r.as_posix().split(':', 1)[1]}"


def kv(out: str) -> dict:
    return {p[0]: p[1] for p in (l.split() for l in out.splitlines()) if len(p) == 2}


def run_ours(tsv: Path) -> dict:
    o = subprocess.run([str(OURS), str(tsv), str(W), str(H), str(FOV)],
                       capture_output=True, text=True).stdout
    return kv(o)


def run_lost(tsv: Path) -> dict:
    cmd = f"{wsl(LOST)} {wsl(tsv)} {wsl(LOSTDB)} {W} {H} {FOV} 0.05 0.0001"
    o = subprocess.run(["wsl", "-e", "bash", "-lc", cmd],
                       capture_output=True, text=True).stdout
    return kv(o)


def score(solved, ra, de, tra, tde):
    if not solved:
        return "REFUSE", None
    e = angular_error_degrees(ra, de, tra, tde)
    return ("TRUE" if e < THR else "FALSE"), e


def main() -> None:
    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    rng = np.random.default_rng(11)
    tally = {m: {s: {"TRUE": 0, "FALSE": 0, "REFUSE": 0} for s in ("ours", "tetra3", "lost")}
             for m in ("clean", "yflip")}
    reasons: dict[str, int] = {}

    print(f"{'case':>4} {'mode':>6} | {'ours':>22} | {'tetra3':>16} | {'LOST':>16}")
    print("-" * 74)

    for i in range(12):
        ra = float(rng.uniform(0, 360))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1, 1))))
        roll = float(rng.uniform(0, 360))
        px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                      fov_deg=FOV, width=W, height=H)
        if len(px) < 8:
            continue
        img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                           noise_sigma=0.5, psf_sigma=1.0, seed=500 + i)
        c = np.asarray(tetra3.get_centroids_from_image(
            Image.fromarray(img), sigma=3, max_returned=30))
        if len(c) < 8:
            continue
        xy = np.column_stack((c[:, 1], c[:, 0]))  # tetra3 gives (y,x)

        for mode in ("clean", "yflip"):
            p = xy.copy()
            if mode == "yflip":
                p[:, 1] = (H - 1) - p[:, 1]
            with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
                tsv = Path(f.name)
                f.write("# x y intensity\n")
                for j, (x, y) in enumerate(p):
                    f.write(f"{x:.6f} {y:.6f} {2000 - j}\n")

            o = run_ours(tsv)
            l = run_lost(tsv)
            t = tetra3.Tetra3().solve_from_centroids(
                np.column_stack((p[:, 1], p[:, 0])), (H, W),
                fov_estimate=FOV, fov_max_error=1.0)
            tsv.unlink(missing_ok=True)

            so, eo = score(o.get("attitude_known") == "1",
                           float(o.get("attitude_ra", 0)), float(o.get("attitude_de", 0)),
                           ra, dec)
            sl, el = score(l.get("attitude_known") == "1",
                           float(l.get("attitude_ra", 0)), float(l.get("attitude_de", 0)),
                           ra, dec)
            st, et = score(t.get("RA") is not None,
                           t.get("RA") or 0.0, t.get("Dec") or 0.0, ra, dec)

            for name, s in (("ours", so), ("tetra3", st), ("lost", sl)):
                tally[mode][name][s] += 1
            r = o.get("gate_reason", "?")
            reasons[r] = reasons.get(r, 0) + 1

            def fmt(s, e):
                return f"{s:>6} {(f'{e:8.3f}' if e is not None else '       -')}"

            print(f"{i:>4} {mode:>6} | {fmt(so, eo)} {o.get('gate_reason','')[:12]:>13} "
                  f"| {fmt(st, et):>16} | {fmt(sl, el):>16}")

    print("-" * 74)
    for mode in ("clean", "yflip"):
        print(f"\n  {mode}:")
        for name in ("ours", "tetra3", "lost"):
            d = tally[mode][name]
            print(f"    {name:>7}:  TRUE {d['TRUE']:>2}   FALSE {d['FALSE']:>2}   REFUSE {d['REFUSE']:>2}")
    print(f"\n  ours gate reasons: {reasons}")


if __name__ == "__main__":
    main()
