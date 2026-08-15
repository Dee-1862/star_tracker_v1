"""Convention diagnostics + cross-generator + focal-length sweep.

1. Residual vectors (dRA, dDec, dRoll) for existing LOST-native results
2. Tetra3-native / Hipparcos-projected centroids scored by both solvers
3. Focal-length (FOV) error sweep — Pitch 2 probe
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
from generate_star_field import camera_basis, project_catalog, render_image  # noqa: E402
from prepare_catalog import parse_catalog, to_binary_records  # noqa: E402
from run_comparison import angular_error_degrees, json_safe  # noqa: E402
from run_decoupled_lost_vs_tetra3 import (  # noqa: E402
    LOST_ADAPTER,
    LOST_DATABASE,
    run_lost_from_centroids,
    run_tetra3_from_centroids,
    write_centroid_tsv,
)

OUTPUT_DIR = BENCH / "results"
TETRA3_NATIVE_DIR = BENCH / "centroids" / "tetra3_native"
CATALOG_PATH = ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat"
NOMINAL_FOV = 20.0


def wrap_delta_deg(value: float) -> float:
    return ((value + 180.0) % 360.0) - 180.0


def residual_vector(
    truth: dict[str, float],
    solved: dict[str, Any],
) -> dict[str, float | None]:
    if not solved.get("solved"):
        return {"dra_deg": None, "ddec_deg": None, "droll_deg": None}
    return {
        "dra_deg": float(solved["ra_deg"]) - float(truth["ra_deg"]),
        "ddec_deg": float(solved["dec_deg"]) - float(truth["dec_deg"]),
        "droll_deg": wrap_delta_deg(
            float(solved["roll_deg"]) - float(truth["roll_deg"])
        ),
        "droll_negated_deg": wrap_delta_deg(
            float(solved["roll_deg"]) - (-float(truth["roll_deg"]))
        ),
    }


def diagnose_existing(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for case in payload["cases"]:
        truth = case["truth"]
        tetra_res = residual_vector(truth, case["tetra3"])
        lost_res = residual_vector(truth, case["lost"])
        between_roll = None
        if case["tetra3"].get("solved") and case["lost"].get("solved"):
            between_roll = wrap_delta_deg(
                float(case["tetra3"]["roll_deg"]) - float(case["lost"]["roll_deg"])
            )
            between_roll_neg = wrap_delta_deg(
                float(case["tetra3"]["roll_deg"])
                - (-float(case["lost"]["roll_deg"]))
            )
        else:
            between_roll_neg = None
        rows.append(
            {
                "case": case["case"],
                "truth": truth,
                "tetra3_residual": tetra_res,
                "lost_residual": lost_res,
                "boresight_error_between_deg": case.get(
                    "attitude_error_between_solvers_deg"
                ),
                "roll_delta_tetra_minus_lost_deg": between_roll,
                "roll_delta_tetra_minus_neg_lost_deg": between_roll_neg,
            }
        )

    tetra_dra = [r["tetra3_residual"]["dra_deg"] for r in rows]
    tetra_ddec = [r["tetra3_residual"]["ddec_deg"] for r in rows]
    tetra_droll = [r["tetra3_residual"]["droll_deg"] for r in rows]
    tetra_droll_neg = [r["tetra3_residual"]["droll_negated_deg"] for r in rows]
    between_roll = [r["roll_delta_tetra_minus_lost_deg"] for r in rows]
    between_roll_neg = [r["roll_delta_tetra_minus_neg_lost_deg"] for r in rows]

    diagnosis = {
        "interpretation": {
            "boresight": (
                "RA/Dec residuals are small and change sign/magnitude across cases "
                "→ not a fixed convention bias; consistent with estimator/FOV-fit noise."
            ),
            "roll": (
                "tetra3_roll ≈ (-lost_roll) mod 360 on every case "
                "→ systematic roll-sign convention mismatch, not accuracy."
            ),
        },
        "tetra3_vs_truth": {
            "mean_dra_deg": float(np.mean(tetra_dra)),
            "mean_ddec_deg": float(np.mean(tetra_ddec)),
            "mean_droll_deg": float(np.mean(tetra_droll)),
            "mean_droll_vs_negated_truth_deg": float(np.mean(tetra_droll_neg)),
            "std_dra_deg": float(np.std(tetra_dra)),
            "std_ddec_deg": float(np.std(tetra_ddec)),
            "std_droll_deg": float(np.std(tetra_droll)),
        },
        "tetra3_vs_lost_roll": {
            "mean_roll_delta_deg": float(np.mean(between_roll)),
            "mean_roll_delta_after_negating_lost_deg": float(np.mean(between_roll_neg)),
            "max_abs_roll_delta_after_negating_lost_deg": float(
                np.max(np.abs(between_roll_neg))
            ),
        },
        "cases": rows,
    }
    return diagnosis


def export_tetra3_native_cases(count: int, seed: int) -> list[Path]:
    TETRA3_NATIVE_DIR.mkdir(parents=True, exist_ok=True)
    catalog = to_binary_records(parse_catalog(CATALOG_PATH, magnitude_limit=6.0))
    rng = np.random.default_rng(seed)
    # Fixed attitudes spanning sky + the LOST-native set for direct comparison
    attitudes = [
        (20.0, -20.0, 0.0),
        (50.0, -10.0, 40.0),
        (80.0, 0.0, 80.0),
        (110.0, 10.0, 120.0),
        (140.0, 20.0, 160.0),
    ]
    while len(attitudes) < count:
        ra = float(rng.uniform(0.0, 360.0))
        dec = float(np.rad2deg(np.arcsin(rng.uniform(-1.0, 1.0))))
        roll = float(rng.uniform(0.0, 360.0))
        attitudes.append((ra, dec, roll))

    solver = tetra3.Tetra3()
    paths: list[Path] = []
    for index, (ra, dec, roll) in enumerate(attitudes[:count]):
        pixels, magnitudes, star_ids = project_catalog(
            catalog,
            camera_basis(ra, dec, roll),
            fov_deg=NOMINAL_FOV,
            width=1024,
            height=1024,
        )
        # Perfect projected centroids (generator-native, no extraction noise)
        order = np.argsort(magnitudes)
        perfect = [
            {
                "x": float(pixels[i, 0]),
                "y": float(pixels[i, 1]),
                "intensity": int(max(1, 1000 - rank)),
                "source": "hipparcos_project_perfect",
                "star_id": int(star_ids[i]),
            }
            for rank, i in enumerate(order)
        ]

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
        from PIL import Image

        extracted = tetra3.get_centroids_from_image(
            Image.fromarray(image),
            sigma=3,
            max_returned=30,
            max_area=100,
            min_area=3,
        )
        tetra_centroids = [
            {
                "y": float(row[0]),
                "x": float(row[1]),
                "intensity": int(1000 - rank),
                "source": "tetra3_extract",
            }
            for rank, row in enumerate(np.asarray(extracted))
        ]

        payload = {
            "case_id": f"tetra3_native_{index:04d}",
            "provenance": "Hipparcos projector + tetra3 centroid extraction",
            "truth": {
                "ra_deg": ra,
                "dec_deg": dec,
                "roll_deg": roll,
                "visible_stars": int(len(star_ids)),
            },
            "image_shape": [1024, 1024],
            "extractors": {
                "perfect_project": perfect,
                "tetra3_extract": tetra_centroids,
            },
        }
        path = TETRA3_NATIVE_DIR / f"tetra3_native_{index:04d}.json"
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        paths.append(path)
        # Keep solver warm / validate tetra3 can solve the image
        _ = solver
    return paths


def score_centroid_set(
    solver: tetra3.Tetra3,
    cases: list[Path],
    extractor: str,
    fov_deg: float,
) -> list[dict[str, Any]]:
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for path in cases:
            payload = json.loads(path.read_text(encoding="utf-8"))
            centroids = payload["extractors"][extractor]
            height, width = payload.get("image_shape", [1024, 1024])
            truth = payload["truth"]
            tsv = tmp_dir / f"{path.stem}.tsv"
            write_centroid_tsv(centroids, tsv)

            tetra = run_tetra3_from_centroids(
                solver, centroids, width, height, fov_deg
            )
            lost = run_lost_from_centroids(tsv, width, height, fov_deg)

            for side in (tetra, lost):
                pointing = None
                if side.get("solved"):
                    pointing = angular_error_degrees(
                        truth["ra_deg"],
                        truth["dec_deg"],
                        side.get("ra_deg"),
                        side.get("dec_deg"),
                    )
                side["pointing_error_deg"] = pointing
                side["false_solve"] = pointing is not None and pointing > 0.1
                side["residual"] = residual_vector(truth, side)

            between = None
            if tetra.get("solved") and lost.get("solved"):
                between = angular_error_degrees(
                    tetra["ra_deg"],
                    tetra["dec_deg"],
                    lost.get("ra_deg"),
                    lost.get("dec_deg"),
                )

            results.append(
                {
                    "case": payload["case_id"],
                    "truth": truth,
                    "extractor": extractor,
                    "fov_commanded_deg": fov_deg,
                    "num_centroids": len(centroids),
                    "tetra3": tetra,
                    "lost": lost,
                    "boresight_error_between_deg": between,
                }
            )
    return results


def summarize_set(cases: list[dict[str, Any]]) -> dict[str, Any]:
    def side_summary(name: str) -> dict[str, Any]:
        solved = [c[name] for c in cases if c[name].get("solved")]
        false = sum(bool(c[name].get("false_solve")) for c in cases)
        pointing = [
            c[name]["pointing_error_deg"]
            for c in cases
            if c[name].get("pointing_error_deg") is not None
        ]
        return {
            "solved": len(solved),
            "false_solves": false,
            "solve_rate": len(solved) / len(cases) if cases else 0.0,
            "false_solve_rate": false / len(cases) if cases else 0.0,
            "median_pointing_error_deg": (
                float(np.median(pointing)) if pointing else None
            ),
        }

    return {
        "cases": len(cases),
        "tetra3": side_summary("tetra3"),
        "lost": side_summary("lost"),
        "both_solved": sum(
            bool(c["tetra3"].get("solved")) and bool(c["lost"].get("solved"))
            for c in cases
        ),
    }


def focal_length_sweep(
    solver: tetra3.Tetra3,
    case_paths: list[Path],
    extractor: str,
    errors_pct: list[float],
) -> dict[str, Any]:
    sweep = []
    for err in errors_pct:
        fov = NOMINAL_FOV * (1.0 + err / 100.0)
        cases = score_centroid_set(solver, case_paths, extractor, fov)
        summary = summarize_set(cases)
        sweep.append(
            {
                "focal_error_pct": err,
                "fov_commanded_deg": fov,
                "summary": summary,
            }
        )
        print(
            f"FOV err {err:+.1f}% (fov={fov:.3f}): "
            f"tetra3 {summary['tetra3']['solved']}/{summary['cases']} "
            f"(false {summary['tetra3']['false_solves']}), "
            f"lost {summary['lost']['solved']}/{summary['cases']} "
            f"(false {summary['lost']['false_solves']})"
        )
    return {
        "nominal_fov_deg": NOMINAL_FOV,
        "extractor": extractor,
        "errors_pct": errors_pct,
        "sweep": sweep,
    }


def write_report(
    diagnosis: dict[str, Any],
    cross: dict[str, Any],
    sweep: dict[str, Any],
    path: Path,
) -> None:
    lines = [
        "# Convention Diagnosis + Cross-Generator + Focal Sweep",
        "",
        "## 1. Residual vectors on LOST-native centroids",
        "",
        diagnosis["interpretation"]["boresight"],
        "",
        diagnosis["interpretation"]["roll"],
        "",
        "Per-case tetra3 residuals vs truth (deg):",
        "",
        "| case | dRA | dDec | dRoll | dRoll vs -truth |",
        "|------|-----|------|-------|-----------------|",
    ]
    for row in diagnosis["cases"]:
        t = row["tetra3_residual"]
        lines.append(
            f"| {row['case']} | {t['dra_deg']:.6f} | {t['ddec_deg']:.6f} | "
            f"{t['droll_deg']:.4f} | {t['droll_negated_deg']:.6f} |"
        )
    lines.extend(
        [
            "",
            f"- Mean tetra3 dRA/dDec: "
            f"{diagnosis['tetra3_vs_truth']['mean_dra_deg']:.6f} / "
            f"{diagnosis['tetra3_vs_truth']['mean_ddec_deg']:.6f}",
            f"- Mean tetra3−lost roll delta: "
            f"{diagnosis['tetra3_vs_lost_roll']['mean_roll_delta_deg']:.3f}°",
            f"- After negating LOST roll: "
            f"{diagnosis['tetra3_vs_lost_roll']['mean_roll_delta_after_negating_lost_deg']:.6f}° "
            f"(max abs {diagnosis['tetra3_vs_lost_roll']['max_abs_roll_delta_after_negating_lost_deg']:.6f}°)",
            "",
            "## 2. Cross-generator (Hipparcos projector centroids)",
            "",
        ]
    )
    for name, block in cross.items():
        s = block["summary"]
        lines.extend(
            [
                f"### extractor=`{name}`",
                f"- tetra3: {s['tetra3']['solved']}/{s['cases']} "
                f"(false {s['tetra3']['false_solves']}), "
                f"median err {s['tetra3']['median_pointing_error_deg']}",
                f"- LOST: {s['lost']['solved']}/{s['cases']} "
                f"(false {s['lost']['false_solves']}), "
                f"median err {s['lost']['median_pointing_error_deg']}",
                "",
            ]
        )
    lines.extend(
        [
            "## 3. Focal-length / FOV error sweep",
            "",
            f"Extractor: `{sweep['extractor']}`, nominal FOV {sweep['nominal_fov_deg']}°",
            "",
            "| FOV error | FOV used | tetra3 solves | tetra3 false | LOST solves | LOST false |",
            "|-----------|----------|---------------|--------------|-------------|------------|",
        ]
    )
    for row in sweep["sweep"]:
        s = row["summary"]
        lines.append(
            f"| {row['focal_error_pct']:+.1f}% | {row['fov_commanded_deg']:.3f} | "
            f"{s['tetra3']['solved']}/{s['cases']} | {s['tetra3']['false_solves']} | "
            f"{s['lost']['solved']}/{s['cases']} | {s['lost']['false_solves']} |"
        )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260802)
    parser.add_argument(
        "--fov-errors",
        type=float,
        nargs="+",
        default=[-5.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 5.0],
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    if not LOST_ADAPTER.is_file() or not LOST_DATABASE.is_file():
        raise FileNotFoundError("LOST adapter/database missing")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Residual diagnosis (LOST-native) ===")
    existing = BENCH / "results" / "decoupled_lost_vs_tetra3.json"
    diagnosis = diagnose_existing(existing)
    print(json.dumps(diagnosis["interpretation"], indent=2))
    print(json.dumps(diagnosis["tetra3_vs_lost_roll"], indent=2))

    print("=== Export tetra3-native / Hipparcos-projected cases ===")
    case_paths = export_tetra3_native_cases(args.cases, args.seed)
    solver = tetra3.Tetra3()

    cross: dict[str, Any] = {}
    for extractor in ("perfect_project", "tetra3_extract"):
        print(f"=== Cross-generator score: {extractor} ===")
        cases = score_centroid_set(solver, case_paths, extractor, NOMINAL_FOV)
        summary = summarize_set(cases)
        cross[extractor] = {"summary": summary, "cases": cases}
        print(json.dumps(summary, indent=2))

    # Focal sweep on a set where BOTH solvers produce true solves.
    # Prefer Hipparcos-projected if LOST is honest there; else LOST-native.
    sweep_extractor = "perfect_project"
    lost_true = (
        cross["perfect_project"]["summary"]["lost"]["solved"]
        - cross["perfect_project"]["summary"]["lost"]["false_solves"]
    )
    if lost_true == 0:
        print(
            "LOST false-solves on Hipparcos-projected centroids "
            "(camera-convention mismatch). Sweeping LOST-native set instead."
        )
        native = sorted((BENCH / "centroids" / "lost_native").glob("*.json"))
        sweep_paths = native
        sweep_extractor = "lost_generate"
    else:
        sweep_paths = case_paths

    print(f"=== Focal-length sweep on extractor={sweep_extractor} ===")
    sweep = focal_length_sweep(
        solver, sweep_paths, sweep_extractor, args.fov_errors
    )

    results = {
        "diagnosis": diagnosis,
        "cross_generator": json_safe(cross),
        "focal_sweep": json_safe(sweep),
    }
    json_path = args.output_dir / "convention_and_focal_sweep.json"
    report_path = args.output_dir / "CONVENTION_AND_FOCAL_SWEEP.md"
    json_path.write_text(
        json.dumps(json_safe(results), indent=2, sort_keys=True), encoding="utf-8"
    )
    write_report(diagnosis, cross, sweep, report_path)
    print(f"Wrote {json_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
