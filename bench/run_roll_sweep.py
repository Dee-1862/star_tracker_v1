"""HEALPix-uniform sky sweep checking BOTH axes: boresight and roll.

Roll was never scored in earlier runs. This settles (a) whether ours recovers
tilt at all, and (b) whether tetra3's apparent success on mirrored input is
real or an artefact of scoring only the boresight.
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
for p in ("benchmarking", "sim_environment/src", "benchmarking/tetra3_baseline", "bench"):
    sys.path.insert(0, str(ROOT / p))
np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from generate_star_field import camera_basis, project_catalog, render_image  # noqa: E402
from healpix_ring import nside2npix, pix2ang_ring  # noqa: E402
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402
from run_comparison import angular_error_degrees  # noqa: E402

EXE = str(ROOT / "build-vs" / "flight_software" / "Release" /
          "star_tracker_centroid_runner.exe")
LOST = ROOT / "bench" / "adapters" / "lost_from_centroids"
LOSTDB = ROOT / "bench" / "data" / "lost_database_hip.dat"
W = H = 1024
FOV = 20.0
NSIDE = 2
ROLLS = (0.0, 72.0, 144.0, 216.0, 288.0)
POS_THR = 0.05


def wsl(p: Path) -> str:
    r = p.resolve()
    return f"/mnt/{r.drive.rstrip(':').lower()}{r.as_posix().split(':', 1)[1]}"


def kv(o: str) -> dict:
    return {p[0]: p[1] for p in (l.split() for l in o.splitlines()) if len(p) == 2}


def wrap180(a: float) -> float:
    while a > 180.0:
        a -= 360.0
    while a < -180.0:
        a += 360.0
    return a


def write(xy: np.ndarray) -> Path:
    f = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    f.write("# x y intensity\n")
    for j, (x, y) in enumerate(xy):
        f.write(f"{x:.6f} {y:.6f} {3000 - j}\n")
    f.close()
    return Path(f.name)


def main() -> None:
    cat = to_binary_records(parse_catalog(
        ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0))
    theta, phi = pix2ang_ring(NSIDE, list(range(nside2npix(NSIDE))))
    centres = [(float(np.degrees(p)), float(90.0 - np.degrees(t)))
               for t, p in zip(theta, phi)]

    acc: dict[tuple[str, str], list] = {}

    def record(mode, solver, pos_err, roll_err, solved):
        acc.setdefault((mode, solver), []).append((pos_err, roll_err, solved))

    n = 0
    for (ra, dec) in centres:
        for roll in ROLLS:
            px, mags, _ = project_catalog(cat, camera_basis(ra, dec, roll),
                                          fov_deg=FOV, width=W, height=H)
            if len(px) < 8:
                continue
            img = render_image(px, mags, width=W, height=H, noise_mean=2.0,
                               noise_sigma=0.5, psf_sigma=1.0, seed=4000 + n)
            c = np.asarray(tetra3.get_centroids_from_image(
                Image.fromarray(img), sigma=3, max_returned=30))
            if len(c) < 8:
                continue
            n += 1
            xy = np.column_stack((c[:, 1], c[:, 0]))

            for mode in ("clean", "yflip"):
                p = xy.copy()
                if mode == "yflip":
                    p[:, 1] = (H - 1) - p[:, 1]
                tsv = write(p)

                o = kv(subprocess.run([EXE, str(tsv), str(W), str(H), str(FOV)],
                                      capture_output=True, text=True).stdout)
                cmd = (f"{wsl(LOST)} {wsl(tsv)} {wsl(LOSTDB)} {W} {H} {FOV} "
                       f"0.05 0.0001")
                l = kv(subprocess.run(["wsl", "-e", "bash", "-lc", cmd],
                                      capture_output=True, text=True).stdout)
                t = tetra3.Tetra3().solve_from_centroids(
                    np.column_stack((p[:, 1], p[:, 0])), (H, W),
                    fov_estimate=FOV, fov_max_error=1.0)
                tsv.unlink(missing_ok=True)

                for solver, solved, sra, sdec, sroll in (
                    ("ours", o.get("attitude_known") == "1",
                     o.get("attitude_ra"), o.get("attitude_de"), o.get("attitude_roll")),
                    ("lost", l.get("attitude_known") == "1",
                     l.get("attitude_ra"), l.get("attitude_de"), l.get("attitude_roll")),
                    ("tetra3", t.get("RA") is not None,
                     t.get("RA"), t.get("Dec"), t.get("Roll")),
                ):
                    if not solved:
                        record(mode, solver, None, None, False)
                        continue
                    pe = angular_error_degrees(float(sra), float(sdec), ra, dec)
                    re_ = abs(wrap180(float(sroll) - roll))
                    record(mode, solver, pe, re_, True)

    print(f"\nHEALPix nside={NSIDE}, {len(ROLLS)} rolls -> {n} fields, "
          f"{n * 2} cases per solver\n")
    print(f"{'mode':>6} {'solver':>7} | {'solved':>6} {'pos ok':>7} {'ROLL ok':>8} "
          f"| {'pos p50':>9} {'roll p50':>9} {'roll max':>9}")
    print("-" * 78)
    for mode in ("clean", "yflip"):
        for solver in ("ours", "tetra3", "lost"):
            v = acc.get((mode, solver), [])
            if not v:
                continue
            sol = [x for x in v if x[2]]
            posok = sum(1 for x in sol if x[0] < POS_THR)
            rollok = sum(1 for x in sol if x[1] < POS_THR)
            pp = f"{np.median([x[0] for x in sol]):.4f}" if sol else "-"
            rp = f"{np.median([x[1] for x in sol]):.4f}" if sol else "-"
            rm = f"{max(x[1] for x in sol):.3f}" if sol else "-"
            print(f"{mode:>6} {solver:>7} | {len(sol):>6} {posok:>7} {rollok:>8} "
                  f"| {pp:>9} {rp:>9} {rm:>9}")
        print()


if __name__ == "__main__":
    main()
