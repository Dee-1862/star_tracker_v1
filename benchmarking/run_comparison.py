"""Reproducible comparison of the C++ flight pipeline and ESA tetra3.

The synthetic cases use the real Hipparcos main catalogue as scene truth.  The
real-image cases are the two FLIR Blackfly frames versioned by ESA tetra3.
Results are written as machine-readable JSON and a concise Markdown report.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import psutil
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
SIM_SRC = ROOT / "sim_environment" / "src"
TETRA3_ROOT = ROOT / "benchmarking" / "tetra3_baseline"
TETRA3_DATA = TETRA3_ROOT / "examples" / "data"
CATALOG_PATH = ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat"
DEFAULT_RUNNER = (
    ROOT
    / "build-vs"
    / "flight_software"
    / "Release"
    / "star_tracker_benchmark_runner.exe"
)
OUTPUT_DIR = ROOT / "benchmarking" / "results"

sys.path.insert(0, str(SIM_SRC))
sys.path.insert(0, str(TETRA3_ROOT))

# ESA tetra3 commit f9fa2eb predates NumPy 2, which removed np.math.
np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from generate_star_field import (  # noqa: E402
    camera_basis,
    project_catalog,
    render_image,
)
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402


@dataclass
class CppRun:
    centroids: int
    cleaned_centroids: int
    reported_matches: int
    indexed_stars: int
    indexed_pairs: int
    index_us: float
    centroid_us: float
    filter_us: float
    match_us: float
    centroid_object_bytes: int
    filter_object_bytes: int
    matcher_object_bytes: int
    catalog_bytes: int
    peak_process_bytes: int
    points: list[dict[str, Any]]


def parse_cpp_output(output: str) -> CppRun:
    summary_names = (
        "centroids",
        "cleaned_centroids",
        "reported_matches",
        "indexed_stars",
        "indexed_pairs",
        "index_us",
        "centroid_us",
        "filter_us",
        "match_us",
        "centroid_object_bytes",
        "filter_object_bytes",
        "matcher_object_bytes",
        "catalog_bytes",
        "peak_process_bytes",
    )
    values: dict[str, Any] | None = None
    points: list[dict[str, Any]] = []
    integer_names = set(summary_names) - {
        "index_us",
        "centroid_us",
        "filter_us",
        "match_us",
    }

    for line in output.splitlines():
        fields = line.split(",")
        if fields[0] == "SUMMARY":
            if len(fields) != len(summary_names) + 1:
                raise ValueError(f"Unexpected C++ summary: {line}")
            values = {
                name: (int(value) if name in integer_names else float(value))
                for name, value in zip(summary_names, fields[1:], strict=True)
            }
        elif fields[0] == "POINT":
            points.append(
                {
                    "index": int(fields[1]),
                    "x": float(fields[2]),
                    "y": float(fields[3]),
                    "intensity": int(fields[4]),
                    "star_id": int(fields[5]),
                    "votes": int(fields[6]),
                }
            )

    if values is None:
        raise ValueError(f"C++ runner returned no summary:\n{output}")
    values["points"] = points
    return CppRun(**values)


def run_cpp(runner: Path, image: np.ndarray, threshold: int) -> CppRun:
    if image.shape != (1024, 1024) or image.dtype != np.uint8:
        raise ValueError("C++ input must be a 1024x1024 uint8 image")
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
        return parse_cpp_output(completed.stdout)
    finally:
        raw_path.unlink(missing_ok=True)


def angular_error_degrees(
    expected_ra: float,
    expected_dec: float,
    actual_ra: float | None,
    actual_dec: float | None,
) -> float | None:
    if actual_ra is None or actual_dec is None:
        return None
    ra1, dec1, ra2, dec2 = np.deg2rad(
        [expected_ra, expected_dec, actual_ra, actual_dec]
    )
    cosine = (
        np.sin(dec1) * np.sin(dec2)
        + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    )
    return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))


def assess_cpp_against_truth(
    run: CppRun,
    truth_pixels: np.ndarray,
    truth_ids: np.ndarray,
    maximum_distance: float = 3.0,
) -> dict[str, int | float | bool]:
    correct = 0
    incorrect = 0
    localized = 0
    localization_errors: list[float] = []

    for point in run.points:
        if len(truth_pixels) == 0:
            if point["star_id"] != 0:
                incorrect += 1
            continue
        distances = np.linalg.norm(
            truth_pixels - np.array([point["x"], point["y"]]), axis=1
        )
        nearest = int(np.argmin(distances))
        distance = float(distances[nearest])
        if distance <= maximum_distance:
            localized += 1
            localization_errors.append(distance)
        if point["star_id"] != 0:
            if distance <= maximum_distance and point["star_id"] == int(truth_ids[nearest]):
                correct += 1
            else:
                incorrect += 1

    reported = correct + incorrect
    return {
        "localized_centroids": localized,
        "correct_ids": correct,
        "incorrect_ids": incorrect,
        "false_id_rate": (incorrect / reported) if reported else 0.0,
        "identification_success": correct >= 3 and incorrect == 0,
        "centroid_rmse_px": (
            float(np.sqrt(np.mean(np.square(localization_errors))))
            if localization_errors
            else 0.0
        ),
    }


def tetra3_array_bytes(solver: tetra3.Tetra3) -> int:
    arrays = (
        solver.star_table,
        solver.pattern_catalog,
        solver.pattern_largest_edge,
        solver.star_catalog_IDs,
    )
    return int(sum(array.nbytes for array in arrays if array is not None))


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def run_synthetic_cases(
    solver: tetra3.Tetra3,
    runner: Path,
    catalog: np.ndarray,
    count: int,
    seed: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    cases: list[dict[str, Any]] = []
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

        cpp = run_cpp(runner, image, threshold=8)
        cpp_assessment = assess_cpp_against_truth(cpp, pixels, star_ids)

        wall_begin = time.perf_counter()
        tetra_solution = solver.solve_from_image(
            Image.fromarray(image),
            fov_estimate=20.0,
            fov_max_error=1.0,
            return_matches=True,
            sigma=3,
            max_returned=30,
            max_area=100,
            min_area=3,
        )
        tetra_wall_ms = (time.perf_counter() - wall_begin) * 1000.0
        pointing_error = angular_error_degrees(
            ra,
            dec,
            tetra_solution.get("RA"),
            tetra_solution.get("Dec"),
        )
        tetra_ids = {
            int(star_id)
            for star_id in (tetra_solution.get("matched_catID") or [])
        }
        truth_id_set = {int(star_id) for star_id in star_ids}
        wrong_tetra_ids = len(tetra_ids - truth_id_set)

        cases.append(
            {
                "case": index,
                "truth": {
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "roll_deg": roll,
                    "visible_stars": len(star_ids),
                },
                "cpp": {
                    **asdict(cpp),
                    **cpp_assessment,
                    "pipeline_ms": (
                        cpp.centroid_us + cpp.filter_us + cpp.match_us
                    )
                    / 1000.0,
                },
                "tetra3": {
                    "solved": tetra_solution.get("RA") is not None,
                    "pointing_error_deg": pointing_error,
                    "matches": tetra_solution.get("Matches"),
                    "incorrect_ids": wrong_tetra_ids,
                    "false_solve": (
                        pointing_error is not None and pointing_error > 0.1
                    ),
                    "reported_solve_ms": tetra_solution.get("T_solve"),
                    "reported_extract_ms": tetra_solution.get("T_extract"),
                    "wall_ms": tetra_wall_ms,
                },
            }
        )
    return cases


def stretch_and_pad_real_image(image: np.ndarray) -> tuple[np.ndarray, int]:
    black = float(np.percentile(image, 50.0))
    white = float(np.percentile(image, 99.99))
    stretched = np.clip((image.astype(np.float64) - black) * 255.0 / (white - black), 0, 255)
    stretched = stretched.astype(np.uint8)
    padded = np.zeros((1024, 1024), dtype=np.uint8)
    y_offset = (1024 - stretched.shape[0]) // 2
    padded[y_offset : y_offset + stretched.shape[0], : stretched.shape[1]] = stretched
    return padded, y_offset


def assess_real_agreement(
    cpp: CppRun,
    tetra_solution: dict[str, Any],
    y_offset: int,
) -> dict[str, int | float]:
    tetra_centroids = tetra_solution.get("matched_centroids") or []
    tetra_ids = tetra_solution.get("matched_catID") or []
    if not tetra_centroids or not tetra_ids:
        return {"agreed_ids": 0, "disagreed_ids": cpp.reported_matches, "agreement_rate": 0.0}

    truth_pixels = np.array(
        [[centroid[1], centroid[0] + y_offset] for centroid in tetra_centroids],
        dtype=np.float64,
    )
    assessment = assess_cpp_against_truth(
        cpp, truth_pixels, np.asarray(tetra_ids), maximum_distance=4.0
    )
    reported = int(assessment["correct_ids"]) + int(assessment["incorrect_ids"])
    return {
        "agreed_ids": int(assessment["correct_ids"]),
        "disagreed_ids": int(assessment["incorrect_ids"]),
        "agreement_rate": (
            float(assessment["correct_ids"]) / reported if reported else 0.0
        ),
    }


def run_real_cases(
    solver: tetra3.Tetra3,
    runner: Path,
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for path in sorted(TETRA3_DATA.glob("*.tiff")):
        original = np.asarray(Image.open(path))
        cpp_image, y_offset = stretch_and_pad_real_image(original)
        cpp = run_cpp(runner, cpp_image, threshold=32)

        wall_begin = time.perf_counter()
        tetra_solution = solver.solve_from_image(
            Image.fromarray(original),
            distortion=[-0.2, 0.1],
            return_matches=True,
        )
        tetra_wall_ms = (time.perf_counter() - wall_begin) * 1000.0
        agreement = assess_real_agreement(cpp, tetra_solution, y_offset)

        cases.append(
            {
                "image": path.name,
                "provenance": "ESA tetra3 FLIR Blackfly real-world test frame",
                "cpp": {
                    **asdict(cpp),
                    **agreement,
                    "pipeline_ms": (
                        cpp.centroid_us + cpp.filter_us + cpp.match_us
                    )
                    / 1000.0,
                    "preprocessing": (
                        "median-to-99.99th-percentile 16-to-8-bit stretch; "
                        "vertically centered zero-padding to 1024x1024"
                    ),
                },
                "tetra3": {
                    "solved": tetra_solution.get("RA") is not None,
                    "ra_deg": tetra_solution.get("RA"),
                    "dec_deg": tetra_solution.get("Dec"),
                    "roll_deg": tetra_solution.get("Roll"),
                    "fov_deg": tetra_solution.get("FOV"),
                    "matches": tetra_solution.get("Matches"),
                    "rmse_arcsec": tetra_solution.get("RMSE"),
                    "reported_solve_ms": tetra_solution.get("T_solve"),
                    "reported_extract_ms": tetra_solution.get("T_extract"),
                    "wall_ms": tetra_wall_ms,
                },
            }
        )
    return cases


def summarize(
    synthetic: list[dict[str, Any]],
    real: list[dict[str, Any]],
    solver: tetra3.Tetra3,
    rss_before_solver: int,
    rss_after_solver: int,
) -> dict[str, Any]:
    cpp_core_bytes = 0
    cpp_peak_bytes = 0
    if synthetic:
        first = synthetic[0]["cpp"]
        cpp_core_bytes = (
            first["centroid_object_bytes"]
            + first["filter_object_bytes"]
            + first["matcher_object_bytes"]
            + first["catalog_bytes"]
        )
        cpp_peak_bytes = max(case["cpp"]["peak_process_bytes"] for case in synthetic)

    cpp_reported = sum(case["cpp"]["reported_matches"] for case in synthetic)
    cpp_correct = sum(case["cpp"]["correct_ids"] for case in synthetic)
    cpp_incorrect = sum(case["cpp"]["incorrect_ids"] for case in synthetic)
    tetra_false = sum(bool(case["tetra3"]["false_solve"]) for case in synthetic)

    return {
        "synthetic_cases": len(synthetic),
        "cpp": {
            "field_successes": sum(
                bool(case["cpp"]["identification_success"]) for case in synthetic
            ),
            "correct_ids": cpp_correct,
            "incorrect_ids": cpp_incorrect,
            "false_id_rate": cpp_incorrect / cpp_reported if cpp_reported else 0.0,
            "median_pipeline_ms": float(
                np.median([case["cpp"]["pipeline_ms"] for case in synthetic])
            ),
            "median_match_ms": float(
                np.median([case["cpp"]["match_us"] / 1000.0 for case in synthetic])
            ),
            "core_static_bytes_including_catalog": cpp_core_bytes,
            "peak_host_process_bytes": cpp_peak_bytes,
        },
        "tetra3": {
            "field_successes": sum(
                bool(case["tetra3"]["solved"])
                and not bool(case["tetra3"]["false_solve"])
                for case in synthetic
            ),
            "false_solves": tetra_false,
            "median_pointing_error_deg": float(
                np.median(
                    [
                        case["tetra3"]["pointing_error_deg"]
                        for case in synthetic
                        if case["tetra3"]["pointing_error_deg"] is not None
                    ]
                )
            ),
            "median_wall_ms": float(
                np.median([case["tetra3"]["wall_ms"] for case in synthetic])
            ),
            "database_array_bytes": tetra3_array_bytes(solver),
            "python_rss_before_solver_bytes": rss_before_solver,
            "python_rss_after_solver_bytes": rss_after_solver,
        },
        "real_images": {
            "count": len(real),
            "tetra3_solved": sum(case["tetra3"]["solved"] for case in real),
            "cpp_tetra3_agreed_ids": sum(case["cpp"]["agreed_ids"] for case in real),
            "cpp_tetra3_disagreed_ids": sum(
                case["cpp"]["disagreed_ids"] for case in real
            ),
        },
    }


def write_report(results: dict[str, Any], path: Path) -> None:
    summary = results["summary"]
    cpp = summary["cpp"]
    tetra = summary["tetra3"]
    real = summary["real_images"]
    synthetic_count = summary["synthetic_cases"]
    report = f"""# Star Tracker Benchmark Report

## Scope and provenance

- C++ candidate: this repository's allocation-free centroiding, three-frame kinematic filter, and pair-distance LIS matcher, Release x64 build.
- Baseline: ESA `tetra3` commit `f9fa2eb9a32a5efc529e2d86f0b59f35b1e9028d`.
- Synthetic validation: {synthetic_count} deterministic 20-degree scenes rendered from the real ESA Hipparcos main catalogue (`hip_main.dat`, magnitude <= 6). These are simulated detector images, not real captures.
- Real-image validation: {real["count"]} official ESA tetra3 FLIR Blackfly test frames. This is a small smoke test, not a statistically representative flight qualification campaign.

## Results

- Synthetic field solves: C++ {cpp["field_successes"]}/{synthetic_count}; tetra3 {tetra["field_successes"]}/{synthetic_count}.
- C++ star IDs: {cpp["correct_ids"]} correct, {cpp["incorrect_ids"]} incorrect; false-ID rate {100.0 * cpp["false_id_rate"]:.2f}%.
- tetra3 false solves: {tetra["false_solves"]}; median pointing error {tetra["median_pointing_error_deg"]:.6f} deg.
- Median synthetic end-to-end host latency: C++ {cpp["median_pipeline_ms"]:.3f} ms (index construction excluded); tetra3 {tetra["median_wall_ms"]:.3f} ms.
- Memory: C++ core objects plus catalogue {cpp["core_static_bytes_including_catalog"] / 1024.0:.1f} KiB; C++ host process peak {cpp["peak_host_process_bytes"] / (1024.0 * 1024.0):.1f} MiB. tetra3 database arrays {tetra["database_array_bytes"] / (1024.0 * 1024.0):.1f} MiB; Python RSS after database load {tetra["python_rss_after_solver_bytes"] / (1024.0 * 1024.0):.1f} MiB. The C++ host-process figure includes runtime/OS overhead and the 1 MiB input buffer; the core figure is the relevant static algorithm budget.
- Official real frames solved by tetra3: {real["tetra3_solved"]}/{real["count"]}. C++ IDs agreeing with tetra3: {real["cpp_tetra3_agreed_ids"]}; disagreeing: {real["cpp_tetra3_disagreed_ids"]}.

## Interpretation

The benchmark is now operational and reproducible, but the two systems are not equivalent. tetra3 returns and verifies a full attitude solution using a four-star pattern database. The current C++ pair-voting matcher returns per-centroid IDs without a geometric attitude verification stage. Any nonzero but incorrect C++ IDs are therefore false identifications, not successful lost-in-space solutions.

Real-frame "agreement" uses tetra3's matched Hipparcos IDs as the reference because the two bundled images do not include independent truth files. It demonstrates execution on real sensor data but must not be presented as independent accuracy validation.

Raw per-case measurements and configuration are in `comparison_results.json`.
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-cases", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260802)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    if not args.runner.is_file():
        raise FileNotFoundError(f"Build the benchmark runner first: {args.runner}")
    if not CATALOG_PATH.is_file():
        raise FileNotFoundError(f"Missing Hipparcos source catalogue: {CATALOG_PATH}")
    real_paths = list(TETRA3_DATA.glob("*.tiff"))
    if len(real_paths) != 2:
        raise FileNotFoundError(
            "Expected the two official tetra3 TIFF files in examples/data"
        )

    process = psutil.Process()
    rss_before_solver = process.memory_info().rss
    solver = tetra3.Tetra3()
    rss_after_solver = process.memory_info().rss

    catalog_frame = parse_catalog(CATALOG_PATH, magnitude_limit=6.0)
    catalog = to_binary_records(catalog_frame)
    synthetic = run_synthetic_cases(
        solver, args.runner, catalog, args.synthetic_cases, args.seed
    )
    real = run_real_cases(solver, args.runner)
    results = {
        "methodology": {
            "seed": args.seed,
            "synthetic_fov_deg": 20.0,
            "synthetic_cpp_threshold": 8,
            "real_cpp_threshold": 32,
            "cpp_index_build_excluded_from_pipeline_latency": True,
            "tetra3_revision": "f9fa2eb9a32a5efc529e2d86f0b59f35b1e9028d",
            "tetra3_database": json_safe(solver.database_properties),
        },
        "summary": summarize(
            synthetic, real, solver, rss_before_solver, rss_after_solver
        ),
        "synthetic": synthetic,
        "real": real,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "comparison_results.json"
    report_path = args.output_dir / "REPORT.md"
    json_path.write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_report(results, report_path)
    print(json.dumps(json_safe(results["summary"]), indent=2))
    print(f"Wrote {json_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
