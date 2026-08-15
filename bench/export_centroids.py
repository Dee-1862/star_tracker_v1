"""Export frozen centroid lists for decoupled star-ID benchmarking.

Each output JSON file contains one test case with centroids from multiple extractors so
adapters can be compared on identical inputs.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SIM_SRC = ROOT / "sim_environment" / "src"
TETRA3_ROOT = ROOT / "benchmarking" / "tetra3_baseline"
TETRA3_DATA = TETRA3_ROOT / "examples" / "data"
DEFAULT_RUNNER = (
    ROOT / "build-vs" / "flight_software" / "Release" / "star_tracker_benchmark_runner.exe"
)
OUTPUT_DIR = Path(__file__).resolve().parent / "centroids"

sys.path.insert(0, str(SIM_SRC))
sys.path.insert(0, str(TETRA3_ROOT))

np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from generate_star_field import (  # noqa: E402
    camera_basis,
    project_catalog,
    render_image,
)
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402


def tetra3_centroids(image: np.ndarray, **kwargs: Any) -> list[dict[str, float]]:
    centroids = tetra3.get_centroids_from_image(Image.fromarray(image), **kwargs)
    return [
        {"y": float(row[0]), "x": float(row[1]), "source": "tetra3"}
        for row in np.asarray(centroids)
    ]


def parse_cpp_points(output: str) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    for line in output.splitlines():
        fields = line.split(",")
        if fields[0] != "POINT":
            continue
        points.append(
            {
                "index": int(fields[1]),
                "x": float(fields[2]),
                "y": float(fields[3]),
                "intensity": int(fields[4]),
                "source": "ours_cpp",
            }
        )
    return points


def ours_centroids(runner: Path, image: np.ndarray, threshold: int) -> list[dict[str, float]]:
    with tempfile.NamedTemporaryFile(suffix=".raw", delete=False) as raw_file:
        raw_path = Path(raw_file.name)
        image.tofile(raw_file)
    try:
        completed = subprocess.run(
            [str(runner), str(raw_path), str(threshold)],
            check=True,
            capture_output=True,
            text=True,
        )
        return parse_cpp_points(completed.stdout)
    finally:
        raw_path.unlink(missing_ok=True)


def export_synthetic(
    count: int,
    seed: int,
    runner: Path,
    output_dir: Path,
) -> list[Path]:
    catalog = to_binary_records(
        parse_catalog(ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat", 6.0)
    )
    rng = np.random.default_rng(seed)
    paths: list[Path] = []
    for index in range(count):
        ra = float(rng.uniform(0.0, 360.0))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1.0, 1.0))))
        roll = float(rng.uniform(0.0, 360.0))
        pixels, magnitudes, star_ids = project_catalog(
            catalog,
            camera_basis(ra, dec, roll),
            fov_deg=20.0,
            width=1024,
            height=1024,
        )
        image = render_image(
            pixels,
            magnitudes,
            width=1024,
            height=1024,
            noise_mean=2.0,
            noise_sigma=0.5,
            psf_sigma=1.0,
            seed=seed + index,
        )
        payload = {
            "case_id": f"synthetic_{index:04d}",
            "provenance": "Hipparcos-rendered synthetic 1024x1024",
            "truth": {
                "ra_deg": ra,
                "dec_deg": dec,
                "roll_deg": roll,
                "visible_stars": int(len(star_ids)),
            },
            "image_shape": [1024, 1024],
            "extractors": {
                "tetra3": tetra3_centroids(image, sigma=3, max_returned=30),
                "ours_cpp": ours_centroids(runner, image, threshold=8),
            },
        }
        path = output_dir / f"synthetic_{index:04d}.json"
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        paths.append(path)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-cases", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260802)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    if not args.runner.is_file():
        raise FileNotFoundError(f"Build benchmark runner first: {args.runner}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = export_synthetic(args.synthetic_cases, args.seed, args.runner, args.output_dir)
    print(f"Wrote {len(paths)} centroid files to {args.output_dir}")


if __name__ == "__main__":
    main()
