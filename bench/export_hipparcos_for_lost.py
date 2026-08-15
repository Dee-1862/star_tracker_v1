"""Convert Hipparcos hip_main.dat into LOST's BSC-like TSV catalog format.

LOST CatalogRead expects lines:
  RAJ2000|DEJ2000|NAME|FLAG|MAG_HIGH.MAG_LOW
where NAME is an integer catalog ID and magnitude is split into integer and
fractional parts (hundredths).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "sim_environment" / "src"))

from prepare_catalog import parse_catalog  # noqa: E402


def format_magnitude(vmag: float) -> str:
    # LOST stores magnitude as integer hundredths: high*100 +/- low
    hundredths = int(round(vmag * 100.0))
    sign = "-" if hundredths < 0 else ""
    hundredths = abs(hundredths)
    high = hundredths // 100
    low = hundredths % 100
    return f"{sign}{high}.{low:02d}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hip",
        type=Path,
        default=ROOT / "sim_environment" / "data" / "raw" / "hip_main.dat",
    )
    parser.add_argument("--magnitude-limit", type=float, default=7.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "bench" / "data" / "hipparcos_for_lost.tsv",
    )
    args = parser.parse_args()

    frame = parse_catalog(args.hip, magnitude_limit=args.magnitude_limit)
    lines: list[str] = []
    for row in frame.itertuples(index=False):
        # Hipparcos IDs are integers; RA/Dec in degrees.
        ra = float(row.ra_deg) if hasattr(row, "ra_deg") else float(row[1])
        dec = float(row.dec_deg) if hasattr(row, "dec_deg") else float(row[2])
        star_id = int(row.star_id) if hasattr(row, "star_id") else int(row[0])
        vmag = float(row.vmag) if hasattr(row, "vmag") else float(row[3])
        if not math.isfinite(ra) or not math.isfinite(dec) or not math.isfinite(vmag):
            continue
        # Normalize RA into [0, 360)
        ra = ra % 360.0
        lines.append(
            f"{ra:010.6f}|{dec:+010.6f}|{star_id}| |{format_magnitude(vmag)}"
        )

    # LOST BscParse asserts result.size() > 9000
    if len(lines) <= 9000:
        raise SystemExit(
            f"Only {len(lines)} stars after mag<={args.magnitude_limit}; "
            "raise magnitude limit so LOST's sanity check passes."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(lines)} stars to {args.output}")


if __name__ == "__main__":
    main()
