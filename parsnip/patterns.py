"""Functions and classes to process string data.

As with any text file format, some string manipulation may be required to process CIF
data. The classes and functions in this module provide simple tools for the manipulation
of string data extracted from CIF files by methods in ``parsnip.parse``.

"""
import re

import numpy as np


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
          String that will replace each whitespace. (Default value = "_"")

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

    def __init__(self, patterns: tuple[tuple[str, str]]):
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
