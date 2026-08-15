"""Decoupled star-ID comparison: tetra3 vs LOST on identical frozen centroids.

Exports tetra3 centroids from bench/centroids/*.json to TSV, runs:
  - tetra3.solve_from_centroids
  - lost_from_centroids (pyramid + DQM)

Both see the same (x,y) list. Centroiding is not part of this comparison.
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

ROOT = Path(__file__).resolve().parents[1]
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "benchmarking"))
sys.path.insert(0, str(ROOT / "sim_environment" / "src"))
sys.path.insert(0, str(ROOT / "benchmarking" / "tetra3_baseline"))

np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from run_comparison import angular_error_degrees, json_safe  # noqa: E402

CENTROID_DIR = BENCH / "centroids"
OUTPUT_DIR = BENCH / "results"
LOST_ADAPTER = BENCH / "adapters" / "lost_from_centroids"
LOST_DATABASE = BENCH / "data" / "lost_database_hip.dat"
SYNTHETIC_FOV_DEG = 20.0


def to_wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    rest = resolved.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


def write_centroid_tsv(centroids: list[dict[str, Any]], path: Path) -> None:
    lines: list[str] = ["# x y intensity"]
    for index, row in enumerate(centroids):
        intensity = int(row.get("intensity", 1000 - index))
        lines.append(f"{float(row['x']):.6f} {float(row['y']):.6f} {intensity}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_lost_adapter_output(stdout: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "solved": False,
        "ra_deg": None,
        "dec_deg": None,
        "roll_deg": None,
        "num_centroids": None,
        "num_star_ids": None,
        "star_id_ms": None,
        "attitude_ms": None,
        "total_ms": None,
    }
    for line in stdout.splitlines():
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        key, value = parts[0], parts[1]
        if key == "attitude_known":
            result["solved"] = value not in ("0", "false", "False")
        elif key == "attitude_ra":
            result["ra_deg"] = float(value)
        elif key == "attitude_de":
            result["dec_deg"] = float(value)
        elif key == "attitude_roll":
            result["roll_deg"] = float(value)
        elif key == "num_centroids":
            result["num_centroids"] = int(value)
        elif key == "num_star_ids":
            result["num_star_ids"] = int(value)
        elif key == "starid_average_ns":
            result["star_id_ms"] = float(value) / 1_000_000.0
        elif key == "attitude_average_ns":
            result["attitude_ms"] = float(value) / 1_000_000.0
        elif key == "total_average_ns":
            result["total_ms"] = float(value) / 1_000_000.0
    return result


def run_lost_from_centroids(
    tsv_path: Path,
    width: int,
    height: int,
    fov_deg: float,
    angular_tolerance_deg: float = 0.05,
    max_mismatch_prob: float = 0.0001,
) -> dict[str, Any]:
    if not LOST_ADAPTER.is_file():
        return {
            "available": False,
            "solved": False,
            "reason": "lost_from_centroids binary missing; build with Makefile.lost_adapter",
        }
    if not LOST_DATABASE.is_file():
        return {
            "available": False,
            "solved": False,
            "reason": "LOST database missing",
        }

    cmd = (
        f"{to_wsl_path(LOST_ADAPTER)} "
        f"{to_wsl_path(tsv_path)} "
        f"{to_wsl_path(LOST_DATABASE)} "
        f"{width} {height} {fov_deg} "
        f"{angular_tolerance_deg} {max_mismatch_prob}"
    )
    completed = subprocess.run(
        ["wsl", "-e", "bash", "-lc", cmd],
        capture_output=True,
        text=True,
        check=False,
    )
    parsed = parse_lost_adapter_output(completed.stdout)
    parsed["available"] = True
    parsed["exit_code"] = completed.returncode
    parsed["stderr"] = completed.stderr
    if completed.returncode != 0 and not parsed["solved"]:
        parsed["crashed"] = True
    else:
        parsed["crashed"] = False
    return parsed


def run_tetra3_from_centroids(
    solver: tetra3.Tetra3,
    centroids: list[dict[str, Any]],
    width: int,
    height: int,
    fov_deg: float,
    fov_max_error: float = 1.0,
) -> dict[str, Any]:
    if not centroids:
        return {"solved": False, "reason": "no centroids"}
    array = np.array([[row["y"], row["x"]] for row in centroids], dtype=np.float64)
    solution = solver.solve_from_centroids(
        array,
        (height, width),
        fov_estimate=fov_deg,
        fov_max_error=fov_max_error,
        return_matches=True,
    )
    return {
        "solved": solution.get("RA") is not None,
        "ra_deg": solution.get("RA"),
        "dec_deg": solution.get("Dec"),
        "roll_deg": solution.get("Roll"),
        "fov_deg": solution.get("FOV"),
        "matches": solution.get("Matches"),
        "rmse_arcsec": solution.get("RMSE"),
        "t_solve_ms": solution.get("T_solve"),
        "matched_cat_ids": solution.get("matched_catID"),
    }


def summarize(cases: list[dict[str, Any]], extractor: str) -> dict[str, Any]:
    tetra_solved = sum(bool(case["tetra3"].get("solved")) for case in cases)
    lost_solved = sum(bool(case["lost"].get("solved")) for case in cases)
    tetra_false = sum(bool(case["tetra3"].get("false_solve")) for case in cases)
    lost_false = sum(bool(case["lost"].get("false_solve")) for case in cases)
    agreements = [
        case["attitude_error_between_solvers_deg"]
        for case in cases
        if case.get("attitude_error_between_solvers_deg") is not None
    ]
    tetra_errors = [
        case["tetra3"]["pointing_error_deg"]
        for case in cases
        if case["tetra3"].get("pointing_error_deg") is not None
    ]
    lost_errors = [
        case["lost"]["pointing_error_deg"]
        for case in cases
        if case["lost"].get("pointing_error_deg") is not None
    ]
    return {
        "cases": len(cases),
        "centroid_extractor": extractor,
        "tetra3_solved": tetra_solved,
        "tetra3_false_solves": tetra_false,
        "lost_solved": lost_solved,
        "lost_false_solves": lost_false,
        "median_tetra3_pointing_error_deg": (
            float(np.median(tetra_errors)) if tetra_errors else None
        ),
        "median_lost_pointing_error_deg": (
            float(np.median(lost_errors)) if lost_errors else None
        ),
        "both_solved": sum(
            bool(case["tetra3"].get("solved")) and bool(case["lost"].get("solved"))
            for case in cases
        ),
        "median_attitude_error_between_solvers_deg": (
            float(np.median(agreements)) if agreements else None
        ),
        "lost_database_bytes": (
            LOST_DATABASE.stat().st_size if LOST_DATABASE.is_file() else 0
        ),
    }


def write_report(summary: dict[str, Any], path: Path) -> None:
    report = f"""# Decoupled Star-ID Benchmark: tetra3 vs LOST

## Setup

- Same frozen centroids for both solvers (`extractor={summary["centroid_extractor"]}`)
- FOV estimate: {SYNTHETIC_FOV_DEG} deg
- Cases: {summary["cases"]}
- LOST path: pyramid + DQM via `lost_from_centroids`
- tetra3 path: `solve_from_centroids`

## Results

| Solver | Solves | False solves | Median pointing error |
|--------|--------|--------------|-----------------------|
| tetra3 | {summary["tetra3_solved"]}/{summary["cases"]} | {summary["tetra3_false_solves"]} | {summary["median_tetra3_pointing_error_deg"]} deg |
| LOST   | {summary["lost_solved"]}/{summary["cases"]} | {summary["lost_false_solves"]} | {summary["median_lost_pointing_error_deg"]} deg |

- Both solved: {summary["both_solved"]}/{summary["cases"]}
- Median attitude disagreement (when both solved): {summary["median_attitude_error_between_solvers_deg"]} deg
- LOST database size: {summary["lost_database_bytes"] / 1024.0:.1f} KiB

Raw data: `decoupled_lost_vs_tetra3.json`
"""
    path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centroid-dir", type=Path, default=CENTROID_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--extractor", default="tetra3")
    parser.add_argument("--fov", type=float, default=SYNTHETIC_FOV_DEG)
    args = parser.parse_args()

    paths = sorted(args.centroid_dir.glob("*.json"))
    if not paths:
        raise FileNotFoundError(
            f"No centroid JSON in {args.centroid_dir}. Run export_centroids.py first."
        )
    if not LOST_ADAPTER.is_file():
        raise FileNotFoundError(
            f"Build adapter first: cd bench/adapters && "
            f"make -f Makefile.lost_adapter (inside WSL)"
        )

    solver = tetra3.Tetra3()
    cases: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for path in paths:
            payload = json.loads(path.read_text(encoding="utf-8"))
            centroids = payload["extractors"][args.extractor]
            height, width = payload.get("image_shape", [1024, 1024])
            truth = payload.get("truth", {})

            tsv_path = tmp_dir / f"{path.stem}.tsv"
            write_centroid_tsv(centroids, tsv_path)

            tetra = run_tetra3_from_centroids(
                solver, centroids, width, height, args.fov
            )
            lost = run_lost_from_centroids(tsv_path, width, height, args.fov)

            for side in (tetra, lost):
                pointing = None
                if side.get("solved") and truth:
                    pointing = angular_error_degrees(
                        truth["ra_deg"],
                        truth["dec_deg"],
                        side.get("ra_deg"),
                        side.get("dec_deg"),
                    )
                side["pointing_error_deg"] = pointing
                side["false_solve"] = pointing is not None and pointing > 0.1

            between = None
            if tetra.get("solved") and lost.get("solved"):
                between = angular_error_degrees(
                    tetra["ra_deg"],
                    tetra["dec_deg"],
                    lost.get("ra_deg"),
                    lost.get("dec_deg"),
                )

            cases.append(
                {
                    "case": payload.get("case_id", path.stem),
                    "truth": truth,
                    "num_centroids": len(centroids),
                    "centroid_extractor": args.extractor,
                    "tetra3": tetra,
                    "lost": lost,
                    "attitude_error_between_solvers_deg": between,
                }
            )
            print(
                f"{path.stem}: tetra3="
                f"{'ok' if tetra.get('solved') else 'fail'}"
                f"({tetra.get('pointing_error_deg')}), "
                f"lost={'ok' if lost.get('solved') else 'fail'}"
                f"({lost.get('pointing_error_deg')}), "
                f"delta={between}"
            )

    summary = summarize(cases, args.extractor)
    results = {
        "mode": "decoupled_star_id_lost_vs_tetra3",
        "summary": summary,
        "cases": cases,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "decoupled_lost_vs_tetra3.json"
    report_path = args.output_dir / "DECOUPLED_LOST_VS_TETRA3.md"
    json_path.write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_report(summary, report_path)
    print(json.dumps(json_safe(summary), indent=2))
    print(f"Wrote {json_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
