"""Quick test: axis flips on Hipparcos-projected centroids → LOST vs tetra3."""

from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "benchmarking"))
sys.path.insert(0, str(ROOT / "sim_environment" / "src"))
sys.path.insert(0, str(ROOT / "benchmarking" / "tetra3_baseline"))

np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from generate_star_field import camera_basis, project_catalog  # noqa: E402
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402
from run_comparison import angular_error_degrees  # noqa: E402
from run_decoupled_lost_vs_tetra3 import (  # noqa: E402
    run_lost_from_centroids,
    run_tetra3_from_centroids,
    write_centroid_tsv,
)

CATALOG_PATH = ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat"
FOV = 20.0
WIDTH = 1024
HEIGHT = 1024

ATTITUDES = [
    (20.0, -20.0, 0.0),
    (50.0, -10.0, 40.0),
    (80.0, 0.0, 80.0),
    (110.0, 10.0, 120.0),
    (140.0, 20.0, 160.0),
]


def transform(pixels: np.ndarray, mode: str) -> np.ndarray:
    out = pixels.copy()
    if mode == "none":
        return out
    if mode == "yflip":
        out[:, 1] = (HEIGHT - 1) - out[:, 1]
        return out
    if mode == "xflip":
        out[:, 0] = (WIDTH - 1) - out[:, 0]
        return out
    if mode == "xyflip":
        out[:, 0] = (WIDTH - 1) - out[:, 0]
        out[:, 1] = (HEIGHT - 1) - out[:, 1]
        return out
    if mode == "swapxy":
        return out[:, ::-1].copy()
    raise ValueError(mode)


def score_mode(solver, catalog, mode: str) -> dict:
    tetra_ok = tetra_false = lost_ok = lost_false = 0
    lost_errs = []
    tetra_errs = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for index, (ra, dec, roll) in enumerate(ATTITUDES):
            pixels, magnitudes, _ = project_catalog(
                catalog,
                camera_basis(ra, dec, roll),
                fov_deg=FOV,
                width=WIDTH,
                height=HEIGHT,
            )
            order = np.argsort(magnitudes)
            pix = transform(pixels[order], mode)
            centroids = [
                {
                    "x": float(pix[i, 0]),
                    "y": float(pix[i, 1]),
                    "intensity": int(max(1, 1000 - i)),
                }
                for i in range(len(pix))
            ]
            tsv = tmp_dir / f"{index}.tsv"
            write_centroid_tsv(centroids, tsv)
            tetra = run_tetra3_from_centroids(
                solver, centroids, WIDTH, HEIGHT, FOV
            )
            lost = run_lost_from_centroids(tsv, WIDTH, HEIGHT, FOV)
            for side, ok_name in ((tetra, "tetra"), (lost, "lost")):
                err = None
                if side.get("solved"):
                    err = angular_error_degrees(
                        ra, dec, side.get("ra_deg"), side.get("dec_deg")
                    )
                false = err is not None and err > 0.1
                true = err is not None and err <= 0.1
                if ok_name == "tetra":
                    tetra_ok += int(true)
                    tetra_false += int(false)
                    if err is not None:
                        tetra_errs.append(err)
                else:
                    lost_ok += int(true)
                    lost_false += int(false)
                    if err is not None:
                        lost_errs.append(err)
    return {
        "mode": mode,
        "tetra3_true": tetra_ok,
        "tetra3_false": tetra_false,
        "lost_true": lost_ok,
        "lost_false": lost_false,
        "tetra3_median_err_deg": float(np.median(tetra_errs)) if tetra_errs else None,
        "lost_median_err_deg": float(np.median(lost_errs)) if lost_errs else None,
    }


def main() -> None:
    catalog = to_binary_records(parse_catalog(CATALOG_PATH, magnitude_limit=6.0))
    solver = tetra3.Tetra3()
    results = []
    for mode in ("none", "yflip", "xflip", "xyflip", "swapxy"):
        print(f"=== mode={mode} ===")
        row = score_mode(solver, catalog, mode)
        results.append(row)
        print(json.dumps(row, indent=2))
    out = BENCH / "results" / "axis_flip_probe.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
