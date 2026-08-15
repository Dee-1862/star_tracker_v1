"""tetra3 adapter: frozen centroids in → attitude + star IDs out."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TETRA3_ROOT = ROOT / "benchmarking" / "tetra3_baseline"
sys.path.insert(0, str(TETRA3_ROOT))

np.math = math  # type: ignore[attr-defined]

import tetra3  # noqa: E402


def solve_from_centroid_file(
    centroid_path: Path,
    extractor: str = "tetra3",
    fov_estimate: float = 20.0,
    fov_max_error: float = 1.0,
) -> dict[str, Any]:
    payload = json.loads(centroid_path.read_text(encoding="utf-8"))
    centroids = payload["extractors"][extractor]
    if not centroids:
        return {"solved": False, "reason": "no centroids"}

    array = np.array([[row["y"], row["x"]] for row in centroids], dtype=np.float64)
    height, width = payload.get("image_shape", [1024, 1024])
    solver = tetra3.Tetra3()
    solution = solver.solve_from_centroids(
        array,
        (height, width),
        fov_estimate=fov_estimate,
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
        "matched_cat_ids": solution.get("matched_catID"),
        "t_solve_ms": solution.get("T_solve"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("centroid_json", type=Path)
    parser.add_argument("--extractor", default="tetra3")
    parser.add_argument("--fov", type=float, default=20.0)
    parser.add_argument("--fov-max-error", type=float, default=1.0)
    args = parser.parse_args()
    result = solve_from_centroid_file(
        args.centroid_json,
        args.extractor,
        args.fov,
        args.fov_max_error,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
