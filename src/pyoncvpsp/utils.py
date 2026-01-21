"""Utilities."""
import numpy as np

__all__ = ['build_mesh']

def build_mesh(zz: float) -> np.ndarray:
    """Construct the radial mesh for a given atomic number zz as in ONCVPSP.

    Args:
        zz (float): Atomic number.

    Returns:
        np.ndarray: Logarithmic radial mesh points.
    """
    amesh: float = 1.006
    al: float = np.log(amesh)
    rr1: float = min( 0.0005 / zz, 0.0005 / 10)
    mmax: int = int(np.log(45.0 / rr1) / al)
    return rr1 * np.exp(al * np.arange(mmax))
