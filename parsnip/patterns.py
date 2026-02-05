# Copyright (c) 2025-2026, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

"""Functions and classes to process string data.

As with any text file format, some string manipulation may be required to process CIF
data. The classes and functions in this module provide simple tools for the manipulation
of string data extracted from CIF files by methods in ``parsnip.parse``.

"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Literal, TypeVar

import numpy as np
from numpy.typing import ArrayLike


def _normalize(string: str | None):
    """Normalize a lookup by stripping spaces and matched quotes."""
    if string is None:
        return ""
    string = string.strip()
    if string[0] == '"' and string[-1] == '"':
        string = string.strip('"')
    elif string[0] == "'" and string[-1] == "'":
        string = string.strip("'")
    return re.sub(r"[\s;]", "", string)


def _normalize_hall(string: str | None):
    """Normalize a lookup by stripping spaces and matched quotes."""
    if string is None:
        return ""
    string = string.strip()
    if string[0] == '"' and string[-1] == '"':
        string = string.strip('"')
    elif string[0] == "'" and string[-1] == "'":
        string = string.strip("'")
    return string


with open(Path(__file__).parent / "symops.json") as f:
    # Process to extract the required data, in the specific format we need
    _full_dict = {
        k: v | {"symops": np.asarray(v["symops"])[:, None]}
        for (k, v) in json.load(f).items()
    }
    SYMOPS_BY_HALL = {_normalize_hall(k): v["symops"] for k, v in _full_dict.items()}
    SYMOPS_BY_HM = {
        _normalize(k): v["symops"]
        for v in _full_dict.values()
        for k in (
            # All four variants are present in COD. We could force the default setting
            # to save a bit of memory, but I'd rather have the accuracy
            v["hermann_mauguin_full"],
            v["hermann_mauguin_full"].split(":")[0],
            v["hermann_mauguin_short"],
            v["hermann_mauguin_short"].split(":")[0],
        )
    }
    SYMOPS_BY_INTL = {
        _normalize(v["table_number"]): v["symops"] for k, v in _full_dict.items()
    }


T = TypeVar("T")

ALLOWED_DELIMITERS = [";\n", "'''", '"""']
"""Delimiters allowed for nonsimple (multi-line) data entries."""


_PROG_STAR = "*+" if sys.version_info >= (3, 11) else "*"
"""Progressively match prefix* if available, else greedily match."""
_PROG_PLUS = "++" if sys.version_info >= (3, 11) else "+"
"""Progressively match prefix+ if available, else greedily match."""

_ANY = r"(?s:.)"
"""Match any character, including newlines."""

_CIF_KEY = r"\S"
"""Match any of the valid characters in a CIF key or loop label.

See Table 1, entry "data-name" of dx.doi.org/10.1107/S1600576715021871
"""

_WHITESPACE = "[\t ]"
"""Officially recognized whitespace characters according to the CIF 1.1 and 2.0 specs.

See section 3.2 of dx.doi.org/10.1107/S1600576715021871 for clarification.
"""

_SAFE_STRING_RE = re.compile(r"(\(\d+\))|[^\d\[\]\,\+\-\/\*\.]")
_SAFE_FRACTN_RE = re.compile(r"([-+]?\d*[/.]?\d+)")


def _contains_wildcard(s: str) -> bool:
    return "?" in s or "*" in s


def _flatten_or_none(ls: list[T]):
    """Return the sole element from a list of l=1, None if l=0, else l."""
    return None if not ls else ls[0] if len(ls) == 1 else ls


def _rational_evaluate_array(arr: str) -> list[list[float]]:
    """Evaluate an array over the ring Q%1."""
    from fractions import Fraction

    one = Fraction(1)
    zero = Fraction(0)

    def _parse_expr(expr: str) -> Fraction:
        """Convert a string into valid Rational numbers and sum."""
        expr = expr.strip().replace("--", "+")
        return sum(Fraction(x) for x in _SAFE_FRACTN_RE.findall(expr)) or zero

    return [
        [
            float(_parse_expr(coord) % one)
            for coord in ls.strip("]").strip("[").split(",")
        ]
        for ls in arr.split("],")
    ]


def _sympy_evaluate_array(arr: str) -> list[list[float]]:
    from sympy import Rational, sympify

    one = Rational(1)
    return [
        [
            float(sympify(coord, rational=True, locals={}) % one)
            for coord in ls.strip("]").strip("[").split(",")
        ]
        for ls in arr.split("],")
    ]


def _safe_eval(
    str_input: str,
    x: int | float,
    y: int | float,
    z: int | float,
    *,
    parse_mode: Literal["python_float", "rational", "sympy"] = "python_float",
) -> list[list[float]]:
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
    # Replace x, y, and z with positional format specifiers and then format in values
    substituted_string = (
        str_input.lower()
        .replace("x", "{0}")
        .replace("y", "{1}")
        .replace("z", "{2}")
        .format(x, y, z)
    )

    # Remove any unexpected characters from the string, including precision specifiers.
    safe_string = _SAFE_STRING_RE.sub("", substituted_string)

    if parse_mode == "rational":
        return _rational_evaluate_array(safe_string)
    if parse_mode == "sympy":
        return _sympy_evaluate_array(safe_string)
    if parse_mode == "python_float":
        return eval(safe_string, {"__builtins__": {}}, {})  # noqa: S307
    raise ValueError(f"Unknown parse mode '{parse_mode}' was provided!")


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


def cast_array_to_float(arr: ArrayLike | None, dtype: type = np.float32):
    """Cast a Numpy array to a dtype, pruning significant digits from numerical values.

    Args:
        arr (np.array[str]): Array of data to convert
        dtype (type, optional):
            dtype to cast array to.
            Default value = ``np.float32``

    Returns
    -------
        np.array[dtype]: Array with new dtype and no significant digit information.
    """
    if arr is None:
        return np.array("nan", dtype=dtype)
    if np.array(arr).shape == (0,):
        return np.array((), dtype=dtype)
    arr = [(el if el is not None else "nan") for el in arr]
    # if any(el is None for el in arr):
    #     raise TypeError("Input array contains `None` and cannot be cast!")
    return np.char.partition(arr, "(")[..., 0].astype(dtype)


def _accumulate_nonsimple_data(data_iter, line: str = ""):
    """Accumulate nonsimmple (multi-line) data entries into a single string."""
    delimiter_count = 0
    while _line_is_continued(data_iter.peek(None)):
        while data_iter.peek(None) and delimiter_count < 2:
            buffer = data_iter.peek().replace(" ", "")
            if buffer[:1] == ";" or any(s in buffer for s in ALLOWED_DELIMITERS):
                delimiter_count += 1
            line += next(data_iter)

        if delimiter_count == 2:
            break  # Exit the context
    return line


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


def _line_is_continued(line: str | None):
    if line is None:
        return False
    line_prefix = line.lstrip()[:3]

    return line_prefix[:1] == ";" or line_prefix == "'''" or line_prefix == '"""'


def _try_cast_to_numeric(s: str):
    """Attempt to cast a string to a number, returning the original string if invalid.

    This method attempts to convert to a float first, followed by an int. Precision
    measurements and indicators of significant digits are stripped.
    """
    parsed = re.match(r"(\d+\.?\d*)", s.strip())
    if parsed is None or re.search(r"[^0-9\.\(\)]", s):
        return s

    if "." in parsed.group(0):
        return float(parsed.group(0))
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
    a2 = np.array([l2 * np.cos(gamma), l2 * np.sin(gamma), 0])
    a3 = np.array([l3 * a3x, l3 * a3y, l3 * a3z])

    return np.array([a1, a2, a3])


def _box_from_lengths_and_angles(l1, l2, l3, alpha, beta, gamma):
    lx = l1
    ly = l2 * np.sin(gamma)

    a3y = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)

    lz = l3 * np.sqrt(1 - np.cos(beta) ** 2 - a3y**2)

    a2x = (l3 * l1 * np.cos(beta)) / lx
    b = l2 * np.cos(gamma)
    c = b * a2x + ly * l3 * a3y

    xy = np.cos(gamma) / np.sin(gamma)
    xz = a2x / lz
    yz = (c - b * a2x) / (ly * lz)

    return tuple(float(x) for x in [lx, ly, lz, xy, xz, yz])


def _lookup_symops(cif) -> np.ndarray | None:
    """Look up the symmetry operations for a space group.

    Note that we choose the default setting (as listed in the International Tables) if
    the provided lookup is ambiguous.

    The space group is extracted from the following keys, in descending priority:
    - _space_group_name_Hall         # Unambiguous but not always present
    - _space_group_name_H-M_alt      # Can include extended (or short) HM symbols
    - _symmetry_space_group_name_H-M # Deprecated, ambiguous setting.
    - _space_group_IT_number         # Ambiguous setting
    - _symmetry_Int_Tables_number    # Deprecated, ambiguous setting
    """
    symops = None
    if (hall := cif["_space_group_name_Hall"]) is not None:
        symops = SYMOPS_BY_HALL.get(_normalize_hall(hall))

    if symops is None and (
        hm := cif["_space_group_name_H-M_alt"] or cif["_symmetry_space_group_name_H-M"]
    ):
        symops = SYMOPS_BY_HM.get(_normalize(hm))

    if symops is None and (
        it := cif["_space_group_IT_number"] or cif["_symmetry_Int_Tables_number"]
    ):
        symops = SYMOPS_BY_INTL.get(_normalize(it))
    return symops
