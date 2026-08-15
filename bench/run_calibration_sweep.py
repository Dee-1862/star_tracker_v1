"""Calibration-robustness sweep: LOST vs tetra3 (TRUE / FALSE / REFUSE).

Design (reviewer-facing):
  - Attitudes: HEALPix RING nside × random roll(s); log galactic |b| and stratify.
  - Focal grid (%): fine near 0 so LOST's edge is resolved.
  - tetra3 fov_max_error is a swept axis (not a constant) — refusal boundary should
    track the search window; time cost should grow with it.
  - LOST: for each focal error, sweep angular-tolerance multipliers and
    max-mismatch-prob; report the *best* config (tuned oracle), not defaults.
  - Outcomes never collapsed: TRUE / FALSE / REFUSE at 0.05° (+ 0.01° / 0.2°).
  - Timing: p50/p95/p99/max per configuration.

Chirality control that makes the axis fix publishable lives in probe_axis_flip.py.

Pilot (default): nside=2 (~48 cells) to confirm tetra3 cliffs track fov_max_error
before scaling to nside=8 (~768).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "benchmarking"))
sys.path.insert(0, str(ROOT / "sim_environment" / "src"))
sys.path.insert(0, str(ROOT / "benchmarking" / "tetra3_baseline"))
sys.path.insert(0, str(BENCH))

np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
import astropy.units as u  # noqa: E402
from generate_star_field import camera_basis, project_catalog  # noqa: E402
from healpix_ring import nside2npix, pix2radec_deg  # noqa: E402
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
CATALOG_PATH = ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat"

NOMINAL_FOV_DEG = 20.0
WIDTH = 1024
HEIGHT = 1024
NOMINAL_ANGULAR_TOL_DEG = 0.05

FOCAL_GRID_PCT = [
    -10.0,
    -7.0,
    -5.0,
    -3.0,
    -2.0,
    -1.5,
    -1.0,
    -0.75,
    -0.5,
    -0.25,
    0.0,
    0.25,
    0.5,
    0.75,
    1.0,
    1.5,
    2.0,
    3.0,
    5.0,
    7.0,
    10.0,
]
TETRA_FOV_MAX_ERROR_DEG = [0.5, 1.0, 2.0, 4.0]
LOST_TOL_MULT = [0.5, 1.0, 2.0, 4.0]
LOST_MISMATCH_PROBS = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
PRIMARY_THRESHOLD_DEG = 0.05
SENSITIVITY_THRESHOLDS_DEG = [0.01, 0.05, 0.2]
GALACTIC_B_BINS = [(0.0, 15.0), (15.0, 30.0), (30.0, 60.0), (60.0, 90.0)]


@dataclass(frozen=True)
class AttitudeSample:
    ipix: int
    roll_index: int
    ra_deg: float
    dec_deg: float
    roll_deg: float
    galactic_b_deg: float
    abs_b_bin: str


def abs_b_bin_label(abs_b: float) -> str:
    for lo, hi in GALACTIC_B_BINS:
        if lo <= abs_b < hi or (hi == 90.0 and abs_b <= hi):
            return f"{lo:.0f}-{hi:.0f}"
    return "out-of-range"


def sample_attitudes(nside: int, rolls: int, seed: int) -> list[AttitudeSample]:
    rng = np.random.default_rng(seed)
    npix = nside2npix(nside)
    samples: list[AttitudeSample] = []
    for ipix in range(npix):
        ra, dec = pix2radec_deg(nside, ipix)
        gal = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs").galactic
        b = float(gal.b.deg)
        for roll_index in range(rolls):
            roll = float(rng.uniform(0.0, 360.0))
            samples.append(
                AttitudeSample(
                    ipix=ipix,
                    roll_index=roll_index,
                    ra_deg=float(ra),
                    dec_deg=float(dec),
                    roll_deg=roll,
                    galactic_b_deg=b,
                    abs_b_bin=abs_b_bin_label(abs(b)),
                )
            )
    return samples


def classify_outcome(
    solved: bool, pointing_error_deg: float | None, threshold_deg: float
) -> str:
    if not solved or pointing_error_deg is None:
        return "REFUSE"
    if pointing_error_deg <= threshold_deg:
        return "TRUE"
    return "FALSE"


def outcome_rank(outcome: str) -> int:
    """Higher is better for oracle selection: TRUE > REFUSE > FALSE."""
    return {"TRUE": 2, "REFUSE": 1, "FALSE": 0}[outcome]


def timing_stats(times_ms: list[float]) -> dict[str, float | None]:
    if not times_ms:
        return {"n": 0, "p50": None, "p95": None, "p99": None, "max": None}
    arr = np.asarray(times_ms, dtype=np.float64)
    return {
        "n": int(arr.size),
        "p50": float(np.percentile(arr, 50)),
        "p95": float(np.percentile(arr, 95)),
        "p99": float(np.percentile(arr, 99)),
        "max": float(np.max(arr)),
    }


def wilson_ci(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 0:
        return (float("nan"), float("nan"))
    phat = successes / n
    denom = 1.0 + z * z / n
    centre = phat + z * z / (2.0 * n)
    margin = z * math.sqrt((phat * (1.0 - phat) + z * z / (4.0 * n)) / n)
    return ((centre - margin) / denom, (centre + margin) / denom)


def perfect_centroids(
    catalog: np.ndarray, attitude: AttitudeSample
) -> list[dict[str, Any]]:
    pixels, magnitudes, star_ids = project_catalog(
        catalog,
        camera_basis(attitude.ra_deg, attitude.dec_deg, attitude.roll_deg),
        fov_deg=NOMINAL_FOV_DEG,
        width=WIDTH,
        height=HEIGHT,
    )
    order = np.argsort(magnitudes)
    return [
        {
            "x": float(pixels[i, 0]),
            "y": float(pixels[i, 1]),
            "intensity": int(max(1, 1000 - rank)),
            "star_id": int(star_ids[i]),
        }
        for rank, i in enumerate(order)
    ]


def score_solver(
    result: dict[str, Any],
    truth: AttitudeSample,
    thresholds: Iterable[float],
) -> dict[str, Any]:
    pointing = None
    if result.get("solved"):
        pointing = angular_error_degrees(
            truth.ra_deg,
            truth.dec_deg,
            float(result["ra_deg"]),
            float(result["dec_deg"]),
        )
    outcomes = {
        f"thr_{t:g}": classify_outcome(bool(result.get("solved")), pointing, t)
        for t in thresholds
    }
    return {
        "solved": bool(result.get("solved")),
        "pointing_error_deg": pointing,
        "outcomes": outcomes,
        "primary": outcomes[f"thr_{PRIMARY_THRESHOLD_DEG:g}"],
        "t_ms": result.get("t_solve_ms")
        if result.get("t_solve_ms") is not None
        else result.get("total_ms"),
    }


def pick_best_lost(
    candidates: list[dict[str, Any]], threshold_key: str
) -> dict[str, Any] | None:
    if not candidates:
        return None

    def key(row: dict[str, Any]) -> tuple:
        outcome = row["score"]["outcomes"][threshold_key]
        err = row["score"]["pointing_error_deg"]
        err_key = err if err is not None else 1e9
        t = row["score"]["t_ms"] if row["score"]["t_ms"] is not None else 1e9
        return (-outcome_rank(outcome), err_key, t)

    return min(candidates, key=key)


def aggregate_counts(outcomes: list[str]) -> dict[str, Any]:
    n = len(outcomes)
    counts = {k: outcomes.count(k) for k in ("TRUE", "FALSE", "REFUSE")}
    true_ci = wilson_ci(counts["TRUE"], n)
    false_ci = wilson_ci(counts["FALSE"], n)
    refuse_ci = wilson_ci(counts["REFUSE"], n)
    return {
        "n": n,
        "TRUE": counts["TRUE"],
        "FALSE": counts["FALSE"],
        "REFUSE": counts["REFUSE"],
        "true_rate": counts["TRUE"] / n if n else None,
        "false_rate": counts["FALSE"] / n if n else None,
        "refuse_rate": counts["REFUSE"] / n if n else None,
        "true_wilson95": true_ci,
        "false_wilson95": false_ci,
        "refuse_wilson95": refuse_ci,
    }


def run_sweep(args: argparse.Namespace) -> dict[str, Any]:
    if not LOST_ADAPTER.is_file() or not LOST_DATABASE.is_file():
        raise FileNotFoundError("LOST adapter/database missing — build in WSL first")

    catalog = to_binary_records(parse_catalog(CATALOG_PATH, magnitude_limit=6.0))
    attitudes = sample_attitudes(args.nside, args.rolls, args.seed)
    if args.max_attitudes is not None:
        attitudes = attitudes[: args.max_attitudes]

    focal_grid = list(args.focal_pct)
    tetra_fov_max = list(args.tetra_fov_max_error)
    lost_tol_mult = list(args.lost_tol_mult) if args.lost_gates else [1.0]
    lost_mismatch = (
        list(args.lost_mismatch_prob) if args.lost_gates else [1e-4]
    )
    thresholds = list(SENSITIVITY_THRESHOLDS_DEG)
    threshold_key = f"thr_{PRIMARY_THRESHOLD_DEG:g}"

    solver = tetra3.Tetra3()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = args.output_dir / args.checkpoint_name
    if checkpoint_path.exists() and not args.resume:
        checkpoint_path.unlink()

    done_keys: set[str] = set()
    records: list[dict[str, Any]] = []
    if args.resume and checkpoint_path.is_file():
        for line in checkpoint_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            records.append(row)
            done_keys.add(row["key"])

    print(
        f"Attitudes={len(attitudes)} nside={args.nside} rolls={args.rolls} "
        f"focal={len(focal_grid)} tetra_fov_max={tetra_fov_max} "
        f"lost_gates={args.lost_gates} "
        f"({len(lost_tol_mult)}×{len(lost_mismatch)} configs)"
    )

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for att_i, attitude in enumerate(attitudes):
            centroids = perfect_centroids(catalog, attitude)
            if len(centroids) < 4:
                print(f"skip ipix={attitude.ipix}: <4 stars")
                continue
            tsv_path = tmp_dir / f"att_{attitude.ipix}_{attitude.roll_index}.tsv"
            write_centroid_tsv(centroids, tsv_path)

            for focal_pct in focal_grid:
                commanded = NOMINAL_FOV_DEG * (1.0 + focal_pct / 100.0)
                base_key = (
                    f"{attitude.ipix}:{attitude.roll_index}:{focal_pct:g}"
                )
                if base_key in done_keys:
                    continue

                tetra_rows: list[dict[str, Any]] = []
                for fov_max in tetra_fov_max:
                    t0 = time.perf_counter()
                    raw = run_tetra3_from_centroids(
                        solver,
                        centroids,
                        WIDTH,
                        HEIGHT,
                        commanded,
                        fov_max_error=fov_max,
                    )
                    wall_ms = (time.perf_counter() - t0) * 1000.0
                    if raw.get("t_solve_ms") is None:
                        raw["t_solve_ms"] = wall_ms
                    score = score_solver(raw, attitude, thresholds)
                    tetra_rows.append(
                        {
                            "fov_max_error_deg": fov_max,
                            "score": score,
                            "fov_reported_deg": raw.get("fov_deg"),
                        }
                    )

                lost_rows: list[dict[str, Any]] = []
                for tol_mult in lost_tol_mult:
                    tol = NOMINAL_ANGULAR_TOL_DEG * tol_mult
                    for mismatch in lost_mismatch:
                        raw = run_lost_from_centroids(
                            tsv_path,
                            WIDTH,
                            HEIGHT,
                            commanded,
                            angular_tolerance_deg=tol,
                            max_mismatch_prob=mismatch,
                        )
                        score = score_solver(raw, attitude, thresholds)
                        lost_rows.append(
                            {
                                "angular_tolerance_deg": tol,
                                "tol_mult": tol_mult,
                                "max_mismatch_prob": mismatch,
                                "score": score,
                            }
                        )

                best_lost = pick_best_lost(lost_rows, threshold_key)
                record = {
                    "key": base_key,
                    "attitude": asdict(attitude),
                    "focal_error_pct": focal_pct,
                    "fov_commanded_deg": commanded,
                    "n_centroids": len(centroids),
                    "tetra3": tetra_rows,
                    "lost": lost_rows,
                    "lost_best": best_lost,
                }
                records.append(record)
                done_keys.add(base_key)
                with checkpoint_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(json_safe(record)) + "\n")

            if (att_i + 1) % max(1, args.progress_every) == 0:
                print(f"  attitudes {att_i + 1}/{len(attitudes)}")

    summary = summarize_records(
        records,
        focal_grid,
        tetra_fov_max,
        thresholds,
        threshold_key,
        args.lost_gates,
    )
    payload = {
        "meta": {
            "nside": args.nside,
            "rolls": args.rolls,
            "seed": args.seed,
            "n_attitudes": len(attitudes),
            "nominal_fov_deg": NOMINAL_FOV_DEG,
            "focal_grid_pct": focal_grid,
            "tetra_fov_max_error_deg": tetra_fov_max,
            "lost_gates": args.lost_gates,
            "lost_tol_mult": lost_tol_mult,
            "lost_mismatch_probs": lost_mismatch,
            "primary_threshold_deg": PRIMARY_THRESHOLD_DEG,
            "sensitivity_thresholds_deg": thresholds,
            "note": (
                "tetra3 fov_max_error is swept; a cliff at ±(fov_max_error/FOV) "
                "is a search-window effect, not an intrinsic algorithm property. "
                "LOST rates use best-of-gates oracle when lost_gates=true."
            ),
        },
        "summary": summary,
    }
    return payload


def summarize_records(
    records: list[dict[str, Any]],
    focal_grid: list[float],
    tetra_fov_max: list[float],
    thresholds: list[float],
    threshold_key: str,
    lost_gates: bool,
) -> dict[str, Any]:
    # tetra3: focal × fov_max_error
    tetra_grid: dict[str, Any] = {}
    for focal in focal_grid:
        tetra_grid[f"{focal:g}"] = {}
        for fov_max in tetra_fov_max:
            outcomes = []
            times = []
            by_b: dict[str, list[str]] = defaultdict(list)
            for rec in records:
                if abs(rec["focal_error_pct"] - focal) > 1e-12:
                    continue
                row = next(
                    r
                    for r in rec["tetra3"]
                    if abs(r["fov_max_error_deg"] - fov_max) < 1e-12
                )
                outcome = row["score"]["outcomes"][threshold_key]
                outcomes.append(outcome)
                if row["score"]["t_ms"] is not None:
                    times.append(float(row["score"]["t_ms"]))
                by_b[rec["attitude"]["abs_b_bin"]].append(outcome)
            tetra_grid[f"{focal:g}"][f"{fov_max:g}"] = {
                "primary": aggregate_counts(outcomes),
                "timing_ms": timing_stats(times),
                "by_abs_b": {
                    bin_name: aggregate_counts(vals) for bin_name, vals in by_b.items()
                },
                "sensitivity": {
                    f"thr_{t:g}": aggregate_counts(
                        [
                            next(
                                r
                                for r in rec["tetra3"]
                                if abs(r["fov_max_error_deg"] - fov_max) < 1e-12
                            )["score"]["outcomes"][f"thr_{t:g}"]
                            for rec in records
                            if abs(rec["focal_error_pct"] - focal) < 1e-12
                        ]
                    )
                    for t in thresholds
                },
            }

    # LOST best-of-gates (or single default)
    lost_best_grid: dict[str, Any] = {}
    for focal in focal_grid:
        outcomes = []
        times = []
        by_b: dict[str, list[str]] = defaultdict(list)
        for rec in records:
            if abs(rec["focal_error_pct"] - focal) > 1e-12:
                continue
            best = rec.get("lost_best")
            if best is None:
                continue
            outcome = best["score"]["outcomes"][threshold_key]
            outcomes.append(outcome)
            if best["score"]["t_ms"] is not None:
                times.append(float(best["score"]["t_ms"]))
            by_b[rec["attitude"]["abs_b_bin"]].append(outcome)
        lost_best_grid[f"{focal:g}"] = {
            "primary": aggregate_counts(outcomes),
            "timing_ms": timing_stats(times),
            "by_abs_b": {
                bin_name: aggregate_counts(vals) for bin_name, vals in by_b.items()
            },
            "oracle": lost_gates,
        }

    # LOST gate response at selected focals (for the prediction test)
    gate_study: dict[str, Any] = {}
    if lost_gates:
        study_focals = [f for f in (-1.0, 0.0, 1.0, 5.0) if f in focal_grid or any(
            abs(f - g) < 1e-12 for g in focal_grid
        )]
        for focal in study_focals:
            # find closest
            focal_use = min(focal_grid, key=lambda g: abs(g - focal))
            by_cfg: dict[str, list[str]] = defaultdict(list)
            for rec in records:
                if abs(rec["focal_error_pct"] - focal_use) > 1e-12:
                    continue
                for row in rec["lost"]:
                    cfg = (
                        f"tol×{row['tol_mult']:g}_mismatch={row['max_mismatch_prob']:g}"
                    )
                    by_cfg[cfg].append(row["score"]["outcomes"][threshold_key])
            gate_study[f"{focal_use:g}"] = {
                cfg: aggregate_counts(vals) for cfg, vals in sorted(by_cfg.items())
            }

    return {
        "tetra3_by_focal_and_fov_max_error": tetra_grid,
        "lost_best_by_focal": lost_best_grid,
        "lost_gate_study": gate_study,
    }


def write_markdown(payload: dict[str, Any], path: Path) -> None:
    meta = payload["meta"]
    summary = payload["summary"]
    lines = [
        "# Calibration robustness sweep",
        "",
        f"- HEALPix nside={meta['nside']}, rolls={meta['rolls']}, "
        f"N≈{meta['n_attitudes']}",
        f"- Primary outcome threshold: {meta['primary_threshold_deg']}° "
        f"(sensitivity {meta['sensitivity_thresholds_deg']})",
        f"- tetra3 `fov_max_error` swept: {meta['tetra_fov_max_error_deg']}",
        f"- LOST oracle gates: {meta['lost_gates']}",
        "",
        meta["note"],
        "",
        "## tetra3 TRUE/FALSE/REFUSE (primary thr) — focal × fov_max_error",
        "",
    ]
    fov_max_list = meta["tetra_fov_max_error_deg"]
    header = "| focal % | " + " | ".join(f"fov_max={v:g}°" for v in fov_max_list) + " |"
    sep = "|---------|" + "|".join(["------"] * len(fov_max_list)) + "|"
    lines.extend([header, sep])
    tetra = summary["tetra3_by_focal_and_fov_max_error"]
    for focal in meta["focal_grid_pct"]:
        cells = []
        for fov_max in fov_max_list:
            cell = tetra[f"{focal:g}"][f"{fov_max:g}"]["primary"]
            cells.append(
                f"T{cell['TRUE']}/F{cell['FALSE']}/R{cell['REFUSE']}"
            )
        lines.append(f"| {focal:+g} | " + " | ".join(cells) + " |")

    lines.extend(
        [
            "",
            "## tetra3 timing (ms) at focal=0 — vs fov_max_error",
            "",
            "| fov_max_error | p50 | p95 | p99 | max |",
            "|---------------|-----|-----|-----|-----|",
        ]
    )
    zero = tetra.get("0") or tetra.get("0.0")
    if zero:
        for fov_max in fov_max_list:
            t = zero[f"{fov_max:g}"]["timing_ms"]
            lines.append(
                f"| {fov_max:g} | {t['p50']} | {t['p95']} | {t['p99']} | {t['max']} |"
            )

    lines.extend(
        [
            "",
            "## LOST best-of-gates (primary thr)",
            "",
            "| focal % | TRUE | FALSE | REFUSE | true_rate |",
            "|---------|------|-------|--------|-----------|",
        ]
    )
    lost = summary["lost_best_by_focal"]
    for focal in meta["focal_grid_pct"]:
        cell = lost[f"{focal:g}"]["primary"]
        lines.append(
            f"| {focal:+g} | {cell['TRUE']} | {cell['FALSE']} | {cell['REFUSE']} | "
            f"{cell['true_rate']} |"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="nside=2, 1 roll — validate cliffs before scaling",
    )
    parser.add_argument("--nside", type=int, default=None)
    parser.add_argument("--rolls", type=int, default=None)
    parser.add_argument("--seed", type=int, default=20260813)
    parser.add_argument("--max-attitudes", type=int, default=None)
    parser.add_argument(
        "--focal-pct",
        type=float,
        nargs="+",
        default=None,
    )
    parser.add_argument(
        "--tetra-fov-max-error",
        type=float,
        nargs="+",
        default=TETRA_FOV_MAX_ERROR_DEG,
    )
    parser.add_argument(
        "--lost-gates",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Sweep LOST tolerance/mismatch and report best-of-gates",
    )
    parser.add_argument(
        "--lost-tol-mult",
        type=float,
        nargs="+",
        default=LOST_TOL_MULT,
    )
    parser.add_argument(
        "--lost-mismatch-prob",
        type=float,
        nargs="+",
        default=LOST_MISMATCH_PROBS,
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--checkpoint-name",
        default="calibration_sweep_checkpoint.jsonl",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--progress-every", type=int, default=4)
    parser.add_argument(
        "--tag",
        default=None,
        help="Output filename tag (default: pilot/full)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.pilot:
        args.nside = args.nside or 2
        args.rolls = args.rolls or 1
        args.focal_pct = args.focal_pct or FOCAL_GRID_PCT
        tag = args.tag or "pilot"
    else:
        args.nside = args.nside or 8
        args.rolls = args.rolls or 1
        args.focal_pct = args.focal_pct or FOCAL_GRID_PCT
        tag = args.tag or "full"

    payload = run_sweep(args)
    json_path = args.output_dir / f"calibration_sweep_{tag}.json"
    md_path = args.output_dir / f"CALIBRATION_SWEEP_{tag.upper()}.md"
    # Drop per-trial bulk from summary JSON (lives in checkpoint jsonl)
    json_path.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_markdown(payload, md_path)
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    print(f"Per-trial checkpoint: {args.output_dir / args.checkpoint_name}")


if __name__ == "__main__":
    main()
