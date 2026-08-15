"""Generate a noisy synthetic star-tracker image from the binary catalog."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import astropy.units as u
import cv2
import numpy as np
from astropy.coordinates import SkyCoord


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
DEFAULT_CATALOG = SIM_ROOT / "data" / "processed" / "hipparcos_vmag6.bin"
DEFAULT_OUTPUT = SIM_ROOT / "output" / "synthetic_star_field.png"


def load_catalog(path: Path) -> np.ndarray:
    """Load and validate the flat, little-endian star catalog."""
    with path.open("rb") as catalog_file:
        header = catalog_file.read(HEADER_SIZE)
        if len(header) != HEADER_SIZE:
            raise ValueError("Catalog header is incomplete")

        magic, record_count, record_size = struct.unpack(HEADER_FORMAT, header)
        if magic != MAGIC:
            raise ValueError(f"Unexpected catalog magic: {magic!r}")
        if record_size != RECORD_DTYPE.itemsize:
            raise ValueError(
                f"Catalog record size is {record_size}, expected {RECORD_DTYPE.itemsize}"
            )

        catalog = np.fromfile(catalog_file, dtype=RECORD_DTYPE, count=record_count)

    if len(catalog) != record_count:
        raise ValueError(f"Catalog contains {len(catalog)} of {record_count} records")
    return catalog


def camera_basis(
    ra_deg: float, dec_deg: float, roll_deg: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return camera right, up, and forward unit vectors in ICRS coordinates."""
    boresight = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    forward = np.asarray(boresight.cartesian.xyz.value, dtype=np.float64)

    ra = boresight.ra.to_value(u.rad)
    dec = boresight.dec.to_value(u.rad)
    east = np.array([-np.sin(ra), np.cos(ra), 0.0], dtype=np.float64)
    north = np.array(
        [-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)],
        dtype=np.float64,
    )

    roll = np.deg2rad(roll_deg)
    right = np.cos(roll) * east + np.sin(roll) * north
    up = -np.sin(roll) * east + np.cos(roll) * north
    return right, up, forward


def project_catalog(
    catalog: np.ndarray,
    basis: tuple[np.ndarray, np.ndarray, np.ndarray],
    fov_deg: float,
    width: int,
    height: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cull stars by dot product and project visible vectors with a pinhole model."""
    if not 0.0 < fov_deg < 180.0:
        raise ValueError("Field of view must be between 0 and 180 degrees")
    if width <= 0 or height <= 0:
        raise ValueError("Image dimensions must be positive")

    right, up, forward = basis
    vectors = np.column_stack((catalog["x"], catalog["y"], catalog["z"])).astype(
        np.float64, copy=False
    )

    forward_component = vectors @ forward
    cone_limit = np.cos(np.deg2rad(fov_deg / 2.0))
    visible = forward_component >= cone_limit

    vectors = vectors[visible]
    depths = forward_component[visible]
    magnitudes = catalog["vmag"][visible]
    star_ids = catalog["star_id"][visible]

    focal_length = (width / 2.0) / np.tan(np.deg2rad(fov_deg / 2.0))
    pixels_x = width / 2.0 + focal_length * (vectors @ right) / depths
    # +up matches LOST/tetra3 camera frames (x=boresight, image y increases with -z).
    # The previous minus produced a mirrored focal plane: pairwise distances still matched
    # catalog patterns (mirror-invariant), so Pyramid reported confident false solves.
    pixels_y = height / 2.0 + focal_length * (vectors @ up) / depths

    on_sensor = (
        (pixels_x >= 0.0)
        & (pixels_x < width)
        & (pixels_y >= 0.0)
        & (pixels_y < height)
    )
    pixels = np.column_stack((pixels_x[on_sensor], pixels_y[on_sensor]))
    return pixels, magnitudes[on_sensor], star_ids[on_sensor]


def render_image(
    pixels: np.ndarray,
    magnitudes: np.ndarray,
    width: int,
    height: int,
    noise_mean: float,
    noise_sigma: float,
    psf_sigma: float,
    seed: int | None,
) -> np.ndarray:
    """Render 5x5 Gaussian PSFs and additive Gaussian detector noise."""
    if noise_sigma < 0.0:
        raise ValueError("Noise sigma cannot be negative")
    if psf_sigma <= 0.0:
        raise ValueError("PSF sigma must be positive")

    rng = np.random.default_rng(seed)
    image = rng.normal(noise_mean, noise_sigma, (height, width)).astype(np.float32)

    offsets = np.arange(-2, 3, dtype=np.int32)
    for (pixel_x, pixel_y), magnitude in zip(pixels, magnitudes, strict=True):
        center_x = int(np.floor(pixel_x))
        center_y = int(np.floor(pixel_y))
        xs = center_x + offsets
        ys = center_y + offsets
        valid_x = (xs >= 0) & (xs < width)
        valid_y = (ys >= 0) & (ys < height)
        if not valid_x.any() or not valid_y.any():
            continue

        patch_x = xs[valid_x]
        patch_y = ys[valid_y]
        dx = patch_x.astype(np.float64) - pixel_x
        dy = patch_y.astype(np.float64) - pixel_y
        gaussian = np.exp(
            -(dy[:, None] ** 2 + dx[None, :] ** 2) / (2.0 * psf_sigma**2)
        )

        # A magnitude-6 star peaks at 20 counts; brighter stars scale physically.
        peak_intensity = min(255.0, 20.0 * 10.0 ** (0.4 * (6.0 - magnitude)))
        image[np.ix_(patch_y, patch_x)] += peak_intensity * gaussian

    return np.clip(image, 0.0, 255.0).astype(np.uint8)


def parse_resolution(value: str) -> tuple[int, int]:
    try:
        width_text, height_text = value.lower().split("x", maxsplit=1)
        width, height = int(width_text), int(height_text)
    except (ValueError, AttributeError) as error:
        raise argparse.ArgumentTypeError("Resolution must use WIDTHxHEIGHT") from error
    if width <= 0 or height <= 0:
        raise argparse.ArgumentTypeError("Resolution dimensions must be positive")
    return width, height


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ra", type=float, required=True, help="Boresight RA in degrees")
    parser.add_argument(
        "--dec", type=float, required=True, help="Boresight declination in degrees"
    )
    parser.add_argument("--roll", type=float, default=0.0, help="Camera roll in degrees")
    parser.add_argument("--fov", type=float, default=20.0, help="Horizontal FOV in degrees")
    parser.add_argument("--resolution", type=parse_resolution, default=(1024, 1024))
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--noise-mean", type=float, default=8.0)
    parser.add_argument("--noise-sigma", type=float, default=2.0)
    parser.add_argument("--psf-sigma", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    width, height = args.resolution
    catalog = load_catalog(args.catalog)
    basis = camera_basis(args.ra, args.dec, args.roll)
    pixels, magnitudes, star_ids = project_catalog(
        catalog, basis, args.fov, width, height
    )
    image = render_image(
        pixels,
        magnitudes,
        width,
        height,
        args.noise_mean,
        args.noise_sigma,
        args.psf_sigma,
        args.seed,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(args.output), image):
        raise OSError(f"OpenCV could not write {args.output}")

    print(f"Projected {len(star_ids):,} stars")
    print(f"Saved {width}x{height} image to {args.output}")


if __name__ == "__main__":
    main()
