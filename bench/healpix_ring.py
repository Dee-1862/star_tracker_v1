"""Pure-Python HEALPix RING pixel centers (no healpy dependency).

Matches healpy.pix2ang(nside, ipix, nest=False) for pixel *centers*.
Formulas follow the public-domain chealpix pix2ang_ring_z_phi.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np

PI = math.pi
TWOPI = 2.0 * math.pi


def nside2npix(nside: int) -> int:
    if nside < 1 or (nside & (nside - 1)) != 0:
        raise ValueError(f"nside must be a power of 2, got {nside}")
    return 12 * nside * nside


def pix2ang_ring(nside: int, ipix: int | Iterable[int]) -> tuple[np.ndarray, np.ndarray]:
    """Return (theta, phi) in radians for RING-ordered pixel centers.

    theta: colatitude [0, pi], phi: longitude [0, 2pi).
    """
    nside = int(nside)
    npix = nside2npix(nside)
    ip = np.atleast_1d(np.asarray(ipix, dtype=np.int64))
    if np.any((ip < 0) | (ip >= npix)):
        raise ValueError("ipix out of range")

    nl4 = 4 * nside
    ncap = 2 * nside * (nside - 1)
    fact1 = 1.5 * nside
    fact2 = 3.0 * nside * nside

    theta = np.empty(ip.shape, dtype=np.float64)
    phi = np.empty(ip.shape, dtype=np.float64)

    # North polar cap
    mask = ip < ncap
    if np.any(mask):
        pix = ip[mask].astype(np.float64)
        iring = np.floor(0.5 * (1.0 + np.sqrt(1.0 + 2.0 * pix))).astype(np.int64)
        iphi = (ip[mask] + 1) - 2 * iring * (iring - 1)
        z = 1.0 - (iring * iring) / fact2
        theta[mask] = np.arccos(z)
        phi[mask] = (iphi - 0.5) * PI / (2.0 * iring)

    # Equatorial belt
    mask = (ip >= ncap) & (ip < (npix - ncap))
    if np.any(mask):
        pix = ip[mask] - ncap
        iring = (pix // nl4) + nside
        iphi = (pix % nl4) + 1
        # 1.0 if (iring+nside) odd, else 0.5
        fodd = np.where(((iring + nside) & 1) != 0, 1.0, 0.5)
        z = (2 * nside - iring) / fact1
        theta[mask] = np.arccos(z)
        phi[mask] = (iphi - fodd) * PI / (2.0 * nside)

    # South polar cap
    mask = ip >= (npix - ncap)
    if np.any(mask):
        ip_s = (npix - ip[mask]).astype(np.float64)
        iring = np.floor(0.5 * (1.0 + np.sqrt(2.0 * ip_s - 1.0))).astype(np.int64)
        iphi = 4 * iring + 1 - (ip_s.astype(np.int64) - 2 * iring * (iring - 1))
        z = -1.0 + (iring * iring) / fact2
        theta[mask] = np.arccos(z)
        phi[mask] = (iphi - 0.5) * PI / (2.0 * iring)

    phi = np.mod(phi, TWOPI)
    if np.isscalar(ipix) or (isinstance(ipix, (int, np.integer))):
        return float(theta.ravel()[0]), float(phi.ravel()[0])
    return theta, phi


def pix2radec_deg(nside: int, ipix: int | Iterable[int]) -> tuple[np.ndarray, np.ndarray]:
    """Pixel centers as (RA, Dec) in degrees."""
    theta, phi = pix2ang_ring(nside, ipix)
    ra = np.rad2deg(np.atleast_1d(phi))
    dec = 90.0 - np.rad2deg(np.atleast_1d(theta))
    if np.isscalar(ipix):
        return float(ra[0]), float(dec[0])
    return ra, dec


def self_check() -> None:
    """Spot-check against published healpy.pix2ang values."""
    t, p = pix2ang_ring(16, 1440)
    assert abs(float(t) - 1.5291175943723188) < 1e-9, (t, p)
    assert abs(float(p) - 0.0) < 1e-9, (t, p)

    t, p = pix2ang_ring(16, [1440, 427, 1520, 0, 3068])
    expect_t = np.array([1.52911759, 0.78550497, 1.57079633, 0.05103658, 3.09055608])
    expect_p = np.array([0.0, 0.78539816, 1.61988371, 0.78539816, 0.78539816])
    assert np.allclose(t, expect_t, atol=1e-8), t
    assert np.allclose(p, expect_p, atol=1e-8), p

    t, p = pix2ang_ring(1, 11)
    assert abs(float(t) - 2.30052398) < 1e-7, t
    assert abs(float(p) - 5.49778714) < 1e-7, p

    assert nside2npix(8) == 768


if __name__ == "__main__":
    self_check()
    print("healpix_ring self_check OK")
