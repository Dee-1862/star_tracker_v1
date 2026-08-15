"""Convert the Hipparcos main catalog into a compact C++-compatible binary.

Binary layout (little-endian):
    Header:  8-byte magic ``STRCAT01``, uint32 record count, uint32 record size
    Record:  uint32 star_id, float32 x, float32 y, float32 z, float32 vmag

The header is 16 bytes and every record is 20 bytes, with no padding.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np
import pandas as pd


MAGIC = b"STRCAT01"
HEADER_FORMAT = "<8sII"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
RECORD_DTYPE = np.dtype(
    [
        ("star_id", "<u4"),
        ("x", "<f4"),
        ("y", "<f4"),
        ("z", "<f4"),
        ("vmag", "<f4"),
    ],
    align=False,
)

SIM_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = SIM_ROOT / "data" / "raw" / "hip_main.dat"
DEFAULT_OUTPUT = SIM_ROOT / "data" / "processed" / "hipparcos_vmag6.bin"


def parse_catalog(path: Path, magnitude_limit: float) -> pd.DataFrame:
    """Read required pipe-delimited columns and reject incomplete records."""
    catalog = pd.read_csv(
        path,
        sep="|",
        header=None,
        usecols=[1, 5, 8, 9],
        dtype=str,
        engine="c",
    )
    catalog.columns = ["star_id", "vmag", "ra_deg", "dec_deg"]

    for column in catalog.columns:
        catalog[column] = pd.to_numeric(catalog[column].str.strip(), errors="coerce")

    catalog = catalog.dropna().loc[lambda frame: frame["vmag"] <= magnitude_limit]
    catalog = catalog.astype({"star_id": np.uint32}).sort_values("star_id")

    if catalog.empty:
        raise ValueError("No valid stars remain after parsing and filtering")
    if catalog["star_id"].duplicated().any():
        raise ValueError("Catalog contains duplicate Hipparcos star IDs")

    return catalog.reset_index(drop=True)


def to_binary_records(catalog: pd.DataFrame) -> np.ndarray:
    """Convert spherical coordinates in degrees to Cartesian unit vectors."""
    ra = np.deg2rad(catalog["ra_deg"].to_numpy(dtype=np.float64))
    dec = np.deg2rad(catalog["dec_deg"].to_numpy(dtype=np.float64))
    cos_dec = np.cos(dec)

    records = np.empty(len(catalog), dtype=RECORD_DTYPE)
    records["star_id"] = catalog["star_id"].to_numpy(dtype=np.uint32)
    records["x"] = cos_dec * np.cos(ra)
    records["y"] = cos_dec * np.sin(ra)
    records["z"] = np.sin(dec)
    records["vmag"] = catalog["vmag"].to_numpy(dtype=np.float32)

    norms = np.sqrt(
        records["x"].astype(np.float64) ** 2
        + records["y"].astype(np.float64) ** 2
        + records["z"].astype(np.float64) ** 2
    )
    if not np.allclose(norms, 1.0, rtol=0.0, atol=2e-7):
        raise ValueError("Coordinate conversion produced non-unit vectors")

    return records


def write_binary(path: Path, records: np.ndarray) -> None:
    """Write the catalog atomically and verify its final byte count."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")

    with temporary_path.open("wb") as output:
        output.write(struct.pack(HEADER_FORMAT, MAGIC, len(records), RECORD_DTYPE.itemsize))
        records.tofile(output)

    expected_size = HEADER_SIZE + len(records) * RECORD_DTYPE.itemsize
    actual_size = temporary_path.stat().st_size
    if actual_size != expected_size:
        temporary_path.unlink(missing_ok=True)
        raise OSError(f"Expected {expected_size} output bytes, wrote {actual_size}")

    temporary_path.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a compact Hipparcos catalog for flight software."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--magnitude-limit", type=float, default=6.0)
    args = parser.parse_args()

    catalog = parse_catalog(args.input, args.magnitude_limit)
    records = to_binary_records(catalog)
    write_binary(args.output, records)

    print(f"Wrote {len(records):,} stars to {args.output}")
    print(f"Binary size: {args.output.stat().st_size:,} bytes")
    print(f"Record size: {RECORD_DTYPE.itemsize} bytes")


if __name__ == "__main__":
    main()
