import numpy as np


def _str2num(val: str):
    """Convert a string value to an integer if possible, or a float otherwise."""
    return float(val) if "." in val else int(val)


def _deg2rad(val: float):
    """Convert a value in degrees to one in radians."""
    return val * np.pi / 180
