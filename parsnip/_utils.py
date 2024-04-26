from __future__ import annotations

import numpy as np


def _str2num(val: str):
    """Convert a string value to an integer if possible, or a float otherwise."""
    return float(val) if "." in val else int(val)


def _deg2rad(val: float | int):
    """Convert a value in degrees to one in radians."""
    return val * np.pi / 180


def _get_distances(positions):
    # Get all indices i!=j
    i_indices, j_indices = np.triu_indices(len(positions), k=1)

    # Compute difference vectors
    r_xyz = positions[i_indices] - positions[j_indices]

    # Compute distances from vectors.
    ij_distances = np.einsum("ij,ij->i", r_xyz, r_xyz, optimize="optimal")

    return ij_distances, i_indices, j_indices
