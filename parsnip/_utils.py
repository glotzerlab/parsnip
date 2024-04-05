import numpy as np


class ParseWarning(Warning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class ParseError(RuntimeError):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def _str2num(val: str):
    """Convert a string value to an integer if possible, or a float otherwise."""
    return float(val) if "." in val else int(val)


def _deg2rad(val: float | int):
    """Convert a value in degrees to one in radians."""
    return val * np.pi / 180
