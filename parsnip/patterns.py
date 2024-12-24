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

from parsnip._errors import ParseWarning

def _safe_eval(str_input: str, x: int | float, y: int | float, z: int | float):
    """Attempt to safely evaluate a string of symmetry equivalent positions.

    Python's ``eval`` is notoriously unsafe. While we could evaluate the entire list at
    once, doing so carries some risk. The typical alternative, ``ast.literal_eval``,
    doesnot work because we need to evaluate mathematical operations.

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

    Returns:
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
        print("(duplicate point, number of occurences)")
        [
            print(pt, count)
            for pt, count in zip(pos[unique_indices], unique_counts)
            if count > 1
        ]

    print()

def cast_array_to_float(arr: np.ndarray, dtype: type = np.float32):
    """Cast a Numpy array to a dtype, pruning significant digits from numerical values.

    Args:
        arr (np.array[str]): Array of data to convert
        dtype (type, optional):
            dtype to cast array to.
            Default value = ``np.float32``

    Returns:
        np.array[float]: Array with new dtype and no significant digit information.
    """
    return np.char.partition(arr, "(")[..., 0].astype(dtype)


def remove_nondelimiting_whitespace(string: str, replacement: str = "_") -> str:
    """Remove nondelimiting whitespaces from a string.

    For the purpose of this function (and CIF files in general), nondelimiting
    whitespaces are those that are enclosed either in single or double quotes.

    Args:
        string (str): Input string to process
        replacement (str):
          String that will replace each nondelimiting whitespace.
          Default value = ``"_"``

    Returns:
        str: String with whitespaces replaced with the replacement character.
    """
    in_quotes = False
    new_str = []
    for char in string:
        if in_quotes and char == " ":
            new_str.append(replacement)
            continue
        else:
            new_str.append(char)

        if char == "'" or char == '"':
            in_quotes = not in_quotes
    return "".join(new_str)


class LineCleaner:
    """Simple object to apply a series of regex patterns to a string.

    To intialize a line cleaner, pass in a tuple of strings of the form
    ``(pattern, replacement)``. Patterns are compiled on initialization to accelerate
    future processing.

    Args:
        patterns (tuple[tuple[str,str]]): Tuple of tuples of strings.
            The first item in each tuple is the pattern to match, and the second item is
            what that pattern will be replaced with.
    """

    def __init__(self, patterns: tuple):
        self.patterns, self.replacements = [], []

        # If we only have a single tuple
        if isinstance(patterns[0], str):
            pattern, replacement = patterns
            self.patterns.append(re.compile(pattern))

            self.replacements.append(replacement)
        else:
            for pattern, replacement in patterns:
                self.patterns.append(re.compile(pattern))

                self.replacements.append(replacement)

    def __call__(self, line: str):
        """Apply patterns defined on initialization of the object to the string.

        ``re.sub(pattern,line)`` is run for each pattern (in order) in self.patterns,
        which is defined on initialization.

        Args:
            line (str): String to apply patterns to.

        Returns:
            str: The substituted lines.
        """
        for pattern, replacement in zip(self.patterns, self.replacements):
            line = pattern.sub(replacement, line)
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
