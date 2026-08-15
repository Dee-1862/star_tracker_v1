"""Three-way benchmark: tetra3, LOST, and our C++ pipeline."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "benchmarking"))
sys.path.insert(0, str(ROOT / "sim_environment" / "src"))

from run_comparison import (  # noqa: E402
    CATALOG_PATH,
    DEFAULT_RUNNER,
    OUTPUT_DIR,
    TETRA3_DATA,
    angular_error_degrees,
    json_safe,
    run_real_cases,
    run_synthetic_cases,
    summarize,
    write_report,
)
from generate_star_field import camera_basis, project_catalog, render_image  # noqa: E402
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402

LOST_BIN = ROOT / "benchmarking" / "lost_baseline" / "lost"
LOST_RUN = BENCH / "adapters" / "lost_run.sh"
LOST_DATABASE = BENCH / "data" / "lost_database_hip.dat"
THREE_WAY_OUTPUT = BENCH / "results"
SYNTHETIC_PNG_DIR = BENCH / "images" / "synthetic"
REAL_PNG_DIR = BENCH / "images" / "real"
REAL_FOV_DEG = 11.4
SYNTHETIC_FOV_DEG = 20.0


def to_wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    rest = resolved.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


def lost_available() -> bool:
    return LOST_BIN.is_file() and LOST_DATABASE.is_file()


def parse_lost_output(stdout: str, stderr: str) -> dict[str, Any]:
    attitude: dict[str, float | None] = {
        "ra_deg": None,
        "dec_deg": None,
        "roll_deg": None,
    }
    timing: dict[str, float | None] = {
        "centroiding_ms": None,
        "star_id_ms": None,
        "attitude_ms": None,
        "total_ms": None,
    }
    attitude_known = False
    combined = stdout + "\n" + stderr
    for line in combined.splitlines():
        stripped = line.strip()
        parts = stripped.split()
        if len(parts) >= 2:
            key, value = parts[0], parts[1]
            if key == "attitude_known":
                attitude_known = value not in ("0", "false", "False")
            elif key == "attitude_ra":
                attitude["ra_deg"] = float(value)
            elif key == "attitude_de":
                attitude["dec_deg"] = float(value)
            elif key == "attitude_roll":
                attitude["roll_deg"] = float(value)
            elif key == "centroiding_average_ns":
                timing["centroiding_ms"] = float(value) / 1_000_000.0
            elif key == "starid_average_ns":
                timing["star_id_ms"] = float(value) / 1_000_000.0
            elif key == "attitude_average_ns":
                timing["attitude_ms"] = float(value) / 1_000_000.0
            elif key == "total_average_ns":
                timing["total_ms"] = float(value) / 1_000_000.0
        if stripped.startswith("RA "):
            attitude["ra_deg"] = float(stripped.split()[1])
        elif stripped.startswith("Dec "):
            attitude["dec_deg"] = float(stripped.split()[1])
        elif stripped.startswith("Roll "):
            attitude["roll_deg"] = float(stripped.split()[1])

    solved = attitude_known or attitude["ra_deg"] is not None
    return {
        "solved": solved,
        **attitude,
        **timing,
        "stdout": stdout,
        "stderr": stderr,
    }


def run_lost_on_png(png_path: Path, fov_deg: float, wall_ms: float | None = None) -> dict[str, Any]:
    if not lost_available():
        return {"available": False, "solved": False, "reason": "LOST binary or database missing"}

    wsl_png = to_wsl_path(png_path)
    wsl_lost = to_wsl_path(LOST_BIN)
    wsl_db = to_wsl_path(LOST_DATABASE)
    cmd = (
        f"{wsl_lost} pipeline "
        f"--png {wsl_png} "
        f"--fov {fov_deg} "
        f"--centroid-algo cog "
        f"--centroid-mag-filter 5 "
        f"--database {wsl_db} "
        f"--star-id-algo py "
        f"--angular-tolerance 0.05 "
        f"--false-stars-estimate 1000 "
        f"--max-mismatch-probability 0.0001 "
        f"--attitude-algo dqm "
        f"--print-attitude - "
        f"--print-speed -"
    )
    begin = time.perf_counter()
    completed = subprocess.run(
        ["wsl", "-e", "bash", "-lc", cmd],
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed_ms = (time.perf_counter() - begin) * 1000.0
    parsed = parse_lost_output(completed.stdout, completed.stderr)
    parsed["available"] = True
    parsed["exit_code"] = completed.returncode
    parsed["wall_ms"] = wall_ms if wall_ms is not None else elapsed_ms
    parsed["crashed"] = completed.returncode not in (0, 11) and not parsed["solved"]
    if completed.returncode == 11 and not parsed["solved"]:
        parsed["no_solve"] = True
    return parsed


def save_grayscale_png(image: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    array = np.asarray(image)
    if array.dtype != np.uint8:
        array = np.clip(array, 0, 255).astype(np.uint8)
    Image.fromarray(array, mode="L").save(path)


def enrich_synthetic_with_lost(
    synthetic: list[dict[str, Any]],
    catalog: np.ndarray,
    seed: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    for case in synthetic:
        index = case["case"]
        truth = case["truth"]
        ra = truth["ra_deg"]
        dec = truth["dec_deg"]
        roll = truth["roll_deg"]
        pixels, magnitudes, _ = project_catalog(
            catalog,
            camera_basis(ra, dec, roll),
            fov_deg=SYNTHETIC_FOV_DEG,
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
        png_path = SYNTHETIC_PNG_DIR / f"synthetic_{index:04d}.png"
        save_grayscale_png(image, png_path)
        lost = run_lost_on_png(png_path, SYNTHETIC_FOV_DEG)
        pointing_error = None
        if lost.get("solved"):
            pointing_error = angular_error_degrees(
                ra, dec, lost.get("ra_deg"), lost.get("dec_deg")
            )
        lost["pointing_error_deg"] = pointing_error
        lost["false_solve"] = (
            pointing_error is not None and pointing_error > 0.1
        )
        case["lost"] = lost
    return synthetic


def enrich_real_with_lost(real: list[dict[str, Any]]) -> list[dict[str, Any]]:
    REAL_PNG_DIR.mkdir(parents=True, exist_ok=True)
    for case in real:
        tiff_path = TETRA3_DATA / case["image"]
        png_path = REAL_PNG_DIR / f"{tiff_path.stem}_gray.png"
        original = np.asarray(Image.open(tiff_path))
        save_grayscale_png(original, png_path)
        lost = run_lost_on_png(png_path, REAL_FOV_DEG)
        tetra = case.get("tetra3", {})
        if lost.get("solved") and tetra.get("solved"):
            lost["attitude_error_vs_tetra3_deg"] = angular_error_degrees(
                tetra["ra_deg"],
                tetra["dec_deg"],
                lost.get("ra_deg"),
                lost.get("dec_deg"),
            )
        else:
            lost["attitude_error_vs_tetra3_deg"] = None
        case["lost"] = lost
    return real


def summarize_lost(
    synthetic: list[dict[str, Any]],
    real: list[dict[str, Any]],
) -> dict[str, Any]:
    synth_lost = [case.get("lost", {}) for case in synthetic]
    real_lost = [case.get("lost", {}) for case in real]
    available = any(item.get("available") for item in synth_lost + real_lost)
    if not available:
        return {"available": False}

    synth_solved = sum(bool(item.get("solved")) for item in synth_lost)
    synth_false = sum(bool(item.get("false_solve")) for item in synth_lost)
    real_solved = sum(bool(item.get("solved")) for item in real_lost)
    real_no_solve = sum(
        bool(item.get("no_solve")) or (item.get("exit_code") == 11 and not item.get("solved"))
        for item in real_lost
    )
    crashes = sum(bool(item.get("crashed")) for item in synth_lost + real_lost)
    walls = [item["wall_ms"] for item in synth_lost if item.get("wall_ms") is not None]
    pointing_errors = [
        item["pointing_error_deg"]
        for item in synth_lost
        if item.get("pointing_error_deg") is not None
    ]
    tetra_agreements = [
        item["attitude_error_vs_tetra3_deg"]
        for item in real_lost
        if item.get("attitude_error_vs_tetra3_deg") is not None
    ]

    return {
        "available": True,
        "database_bytes": LOST_DATABASE.stat().st_size if LOST_DATABASE.is_file() else 0,
        "synthetic_solved": synth_solved,
        "synthetic_false_solves": synth_false,
        "synthetic_crashes": sum(bool(item.get("crashed")) for item in synth_lost),
        "real_solved": real_solved,
        "real_no_solve": real_no_solve,
        "real_crashes": sum(bool(item.get("crashed")) for item in real_lost),
        "total_crashes": crashes,
        "median_wall_ms": float(np.median(walls)) if walls else None,
        "median_pointing_error_deg": (
            float(np.median(pointing_errors)) if pointing_errors else None
        ),
        "real_median_agreement_vs_tetra3_deg": (
            float(np.median(tetra_agreements)) if tetra_agreements else None
        ),
    }


def write_three_way_report(results: dict[str, Any], path: Path) -> None:
    summary = results["summary"]
    lost = results.get("lost_summary", {})
    cpp = summary["cpp"]
    tetra = summary["tetra3"]
    real = summary["real_images"]
    synthetic_count = summary["synthetic_cases"]

    lost_block = "LOST was not run."
    if lost.get("available"):
        lost_block = f"""- LOST synthetic solves: {lost.get("synthetic_solved", 0)}/{synthetic_count}; false solves: {lost.get("synthetic_false_solves", 0)}; crashes: {lost.get("synthetic_crashes", 0)}.
- LOST real-image solves: {lost.get("real_solved", 0)}/{real["count"]}; median agreement vs tetra3: {lost.get("real_median_agreement_vs_tetra3_deg")} deg.
- LOST median wall time (synthetic): {lost.get("median_wall_ms")} ms; database size: {lost.get("database_bytes", 0) / 1024.0:.1f} KiB."""

    report = f"""# Three-Way Star Tracker Benchmark Report

## Solvers

- **Ours:** allocation-free C++ pair-voting matcher (Release x64)
- **tetra3:** ESA 4-star pattern hash + binomial verify
- **LOST:** pyramid + k-vector + DQM attitude (UWCubeSat)

## Results

### Synthetic ({synthetic_count} cases, 20 deg FOV, Hipparcos truth)

| Solver | Field solves | False solves | Median latency |
|--------|--------------|--------------|----------------|
| Ours   | {cpp["field_successes"]}/{synthetic_count} | n/a (86.75% false IDs) | {cpp["median_pipeline_ms"]:.3f} ms |
| tetra3 | {tetra["field_successes"]}/{synthetic_count} | {tetra["false_solves"]} | {tetra["median_wall_ms"]:.3f} ms |
| LOST   | {lost.get("synthetic_solved", "n/a")}/{synthetic_count} | {lost.get("synthetic_false_solves", "n/a")} | {lost.get("median_wall_ms", "n/a")} ms |

### Real images ({real["count"]} ESA FLIR frames)

- tetra3 solved: {real["tetra3_solved"]}/{real["count"]}
- Ours IDs agreeing with tetra3: {real["cpp_tetra3_agreed_ids"]} (disagreed: {real["cpp_tetra3_disagreed_ids"]})
- LOST: see lost_summary in JSON

## LOST detail

{lost_block}

Raw per-case data: `three_way_results.json`
"""
    path.write_text(report, encoding="utf-8")


def run_decoupled_tetra3(centroid_dir: Path) -> list[dict[str, Any]]:
    adapter = BENCH / "adapters" / "tetra3_from_centroids.py"
    results: list[dict[str, Any]] = []
    for path in sorted(centroid_dir.glob("*.json")):
        completed = subprocess.run(
            [sys.executable, str(adapter), str(path), "--extractor", "tetra3"],
            capture_output=True,
            text=True,
            check=False,
        )
        payload = json.loads(path.read_text(encoding="utf-8"))
        solution = json.loads(completed.stdout) if completed.stdout.strip() else {}
        truth = payload.get("truth", {})
        pointing_error = None
        if solution.get("solved") and truth:
            pointing_error = angular_error_degrees(
                truth["ra_deg"],
                truth["dec_deg"],
                solution.get("ra_deg"),
                solution.get("dec_deg"),
            )
        results.append(
            {
                "case": payload.get("case_id", path.stem),
                "tetra3": {**solution, "pointing_error_deg": pointing_error},
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-cases", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260802)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-dir", type=Path, default=THREE_WAY_OUTPUT)
    parser.add_argument("--decoupled", action="store_true")
    parser.add_argument("--skip-lost", action="store_true")
    args = parser.parse_args()

    if args.decoupled:
        centroid_dir = BENCH / "centroids"
        if not any(centroid_dir.glob("*.json")):
            raise FileNotFoundError(
                f"No centroid files in {centroid_dir}. Run bench/export_centroids.py first."
            )
        results = {"mode": "decoupled_star_id", "cases": run_decoupled_tetra3(centroid_dir)}
        args.output_dir.mkdir(parents=True, exist_ok=True)
        out = args.output_dir / "decoupled_results.json"
        out.write_text(json.dumps(json_safe(results), indent=2), encoding="utf-8")
        print(f"Wrote {out}")
        return

    if not args.runner.is_file():
        raise FileNotFoundError(f"Build benchmark runner first: {args.runner}")
    if not CATALOG_PATH.is_file():
        raise FileNotFoundError(f"Missing Hipparcos catalogue: {CATALOG_PATH}")

    import psutil
    import tetra3

    np.math = math  # type: ignore[attr-defined]

    process = psutil.Process()
    rss_before = process.memory_info().rss
    solver = tetra3.Tetra3()
    rss_after = process.memory_info().rss

    catalog = to_binary_records(parse_catalog(CATALOG_PATH, magnitude_limit=6.0))
    synthetic = run_synthetic_cases(
        solver, args.runner, catalog, args.synthetic_cases, args.seed
    )
    real = run_real_cases(solver, args.runner)
    summary = summarize(synthetic, real, solver, rss_before, rss_after)

    lost_summary = {"available": False}
    if not args.skip_lost and lost_available():
        print("Running LOST on synthetic cases...")
        synthetic = enrich_synthetic_with_lost(synthetic, catalog, args.seed)
        print("Running LOST on real images...")
        real = enrich_real_with_lost(real)
        lost_summary = summarize_lost(synthetic, real)
    elif not args.skip_lost:
        lost_summary = {
            "available": False,
            "reason": "LOST binary or database missing",
        }

    results = {
        "mode": "full_pipeline",
        "summary": summary,
        "lost_summary": lost_summary,
        "synthetic": synthetic,
        "real": real,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "three_way_results.json"
    report_path = args.output_dir / "THREE_WAY_REPORT.md"
    json_path.write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_report(results, OUTPUT_DIR / "REPORT.md")
    write_three_way_report(results, report_path)

    print(json.dumps(json_safe({"summary": summary, "lost_summary": lost_summary}), indent=2))
    print(f"Wrote {json_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
