# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Functions and classes to process string data.

As with any text file format, some string manipulation may be required to process CIF
data. The classes and functions in this module provide simple tools for the manipulation
of string data extracted from CIF files by methods in ``parsnip.parse``.

"""
from __future__ import annotations

import re
import warnings

import numpy as np
from numpy.typing import ArrayLike

from parsnip._errors import ParseWarning


def _safe_eval(str_input: str, x: int | float, y: int | float, z: int | float):
    """Attempt to safely evaluate a string of symmetry equivalent positions.

    Python's ``eval`` is notoriously unsafe. While we could evaluate the entire list at
    once, doing so carries some risk. The typical alternative, ``ast.literal_eval``,
    does not work because we need to evaluate mathematical operations.

    We first replace the x,y,z values with ordered fstring inputs, to simplify the input
    of fractional coordinate data. This is done for convenience more than security.

    Once we substitute in the x,y,z values, we should have a string version of a list
    containing only numerics and math operators. We apply a substitution to ensure this
    is the case, then perform one final check. If it passes, we evaluate the list. Note
    that __builtins__ is set to {}, meaning importing functions is not possible. The
    __locals__ dict is also set to {}, so no variables are accessible in the evaluation.

    I cannot guarantee this is fully safe, but it at the very least makes it extremely
    difficult to do any funny business.

    Args:
        str_input (str): String to be evaluated.
        x (int|float): Fractional coordinate in :math:`x`.
        y (int|float): Fractional coordinate in :math:`y`.
        z (int|float): Fractional coordinate in :math:`z`.

    Returns
    -------
        list[list[int|float,int|float,int|float]]:
            :math:`(N,3)` list of fractional coordinates.

    """
    ordered_inputs = {"x": "{0:.20f}", "y": "{1:.20f}", "z": "{2:.20f}"}
    # Replace any x, y, or z with the same character surrounded by curly braces. Then,
    # perform substitutions to insert the actual values.
    substituted_string = (
        re.sub(r"([xyz])", r"{\1}", str_input).format(**ordered_inputs).format(x, y, z)
    )

    # Remove any unexpected characters from the string.
    safe_string = re.sub(r"[^\d\[\]\,\+\-\/\*\.]", "", substituted_string)
    # Double check to be sure:
    assert all(char in ",.0123456789+-/*[]" for char in safe_string), (
        "Evaluation aborted. Check that symmetry operation string only contains "
        "numerics or characters in { [],.+-/ } and adjust `regex_filter` param "
        "accordingly."
    )
    return eval(safe_string, {"__builtins__": {}}, {})  # noqa: S307


def _write_debug_output(unique_indices, unique_counts, pos, check="Initial"):
    print(f"{check} uniqueness check:")
    if len(unique_indices) == len(pos):
        print("... all points are unique (within tolerance).")
    else:
        print("(duplicate point, number of occurrences)")
        [
            print(pt, count)
            for pt, count in zip(np.asarray(pos)[unique_indices], unique_counts)
            if count > 1
        ]

    print()


def cast_array_to_float(arr: ArrayLike, dtype: type = np.float32):
    """Cast a Numpy array to a dtype, pruning significant digits from numerical values.

    Args:
        arr (np.array[str]): Array of data to convert
        dtype (type, optional):
            dtype to cast array to.
            Default value = ``np.float32``

    Returns
    -------
        np.array[float]: Array with new dtype and no significant digit information.
    """
    return np.char.partition(arr, "(")[..., 0].astype(dtype)


def _accumulate_nonsimple_data(data_iter, line=""):
    """Accumulate nonsimmple (multi-line) data entries into a single string."""
    while _line_is_continued(data_iter.peek(None)):
        line += _strip_comments(next(data_iter))
    return _semicolon_to_string(line)


def _is_key(line: str | None):
    return line is not None and line.strip()[:1] == "_"


def _is_data(line: str | None):
    return line is not None and line.strip()[:1] != "_" and line.strip()[:5] != "loop_"


def _strip_comments(s: str):
    return s.split("#")[0].strip()


def _strip_quotes(s: str):
    return s.replace("'", "").replace('"', "")


def _dtype_from_int(i: int):
    return f"<U{i}"


def _semicolon_to_string(line: str):
    if "'" in line and '"' in line:
        warnings.warn(
            (
                "String contains single and double quotes - "
                "line may be parsed incorrectly"
            ),
            ParseWarning,
            stacklevel=2,
        )
    # WARNING: because we split our string, we strip "\n" implicitly
    # This is technically against spec, but is almost never meaningful
    return line.replace(";", "'" if "'" not in line else '"')


def _line_is_continued(line: str | None):
    return line is not None and line.strip()[:1] == ";"


def _try_cast_to_numeric(s: str):
    """Attempt to cast a string to a number, returning the original string if invalid.

    This method attempts to convert to a float first, followed by an int. Precision
    measurements and indicators of significant digits are stripped.
    """
    parsed = re.match(r"(\d+\.?\d*)", s.strip())
    if parsed is None or re.search(r"[^0-9\.\(\)]", s):
        return s
    elif "." in parsed.group(0):
        return float(parsed.group(0))
    else:
        return int(parsed.group(0))


def _matrix_from_lengths_and_angles(l1, l2, l3, alpha, beta, gamma):
    a1 = np.array([l1, 0, 0])
    a2 = np.array([l2 * np.cos(gamma), l2 * np.sin(gamma), 0])
    a3x = np.cos(beta)
    a3y = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    under_sqrt = 1 - a3x**2 - a3y**2
    if under_sqrt < 0:
        raise ValueError("The provided angles can not form a valid box.")
    a3z = np.sqrt(under_sqrt)
    a3 = np.array([l3 * a3x, l3 * a3y, l3 * a3z])
    return np.array([a1, a2, a3])
