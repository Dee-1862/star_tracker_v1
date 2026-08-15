"""Convert LOST print-input-centroids dumps into TSV + JSON for decoupled bench."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "centroids" / "lost_native"


def parse_centroids(path: Path) -> list[dict]:
    xs: dict[int, float] = {}
    ys: dict[int, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        match = re.match(r"input_centroid_(\d+)_x\s+(\S+)", stripped)
        if match:
            xs[int(match.group(1))] = float(match.group(2))
            continue
        match = re.match(r"input_centroid_(\d+)_y\s+(\S+)", stripped)
        if match:
            ys[int(match.group(1))] = float(match.group(2))
    centroids = []
    for rank, index in enumerate(sorted(set(xs) & set(ys))):
        intensity = 1000 - rank
        centroids.append(
            {
                "x": xs[index],
                "y": ys[index],
                "intensity": intensity,
                "source": "lost_generate",
            }
        )
    return centroids


def parse_truth(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    truth: dict[str, float] = {}
    for key, dest in (
        ("expected_attitude_ra", "ra_deg"),
        ("expected_attitude_de", "dec_deg"),
        ("expected_attitude_roll", "roll_deg"),
    ):
        match = re.search(rf"{key}\s+(\S+)", text)
        if match:
            truth[dest] = float(match.group(1))
    return truth


def main() -> None:
    for att in sorted(SRC.glob("*.attitude.txt")):
        case = att.name.replace(".attitude.txt", "")
        cent = SRC / f"{case}.centroids.txt"
        if not cent.is_file():
            continue
        centroids = parse_centroids(cent)
        truth = parse_truth(att)
        tsv = SRC / f"{case}.tsv"
        js = SRC / f"{case}.json"
        lines = ["# x y intensity"]
        for row in centroids:
            lines.append(f"{row['x']:.6f} {row['y']:.6f} {row['intensity']}")
        tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
        payload = {
            "case_id": case,
            "provenance": "LOST --generate --generate-centroids-only (Hipparcos catalog)",
            "truth": truth,
            "image_shape": [1024, 1024],
            "extractors": {"lost_generate": centroids},
        }
        js.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"{case}: {len(centroids)} centroids truth={truth}")


if __name__ == "__main__":
    main()
