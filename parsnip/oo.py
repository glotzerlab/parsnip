# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

from __future__ import annotations

from collections.abc import Iterable
import re
import warnings

import numpy as np
from more_itertools import flatten, peekable

from parsnip._errors import ParseWarning
from parsnip.parse import _parsed_line_generator

NONTABLE_LINE_PREFIXES = ("_", "#")


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


class CifFile:
    """Parsed CIF file."""

    def __init__(self, fn: str, cast_values: bool = False):
        """Create a CifFile object from a filename.

        On construction, the entire file is parsed into key-value pairs and data tables.
        Comment lines are ignored.

        """
        self._fn = fn
        self._pairs = {}
        self._tables = []
        # self._table_keys = []
        # self._table_data = []

        self._cpat = {k: re.compile(pattern) for (k, pattern) in self.PATTERNS.items()}
        self._cast_values = cast_values

        with open(fn) as file:
            self._data = file.read()
        self._parse()


    @property
    def cast_values(self):
        """Whether to cast "number-like" values to ints & floats.

        When set to `True` after construction, the values are modified in-place. This
        action cannot be reversed
        """
        return self._cast_values

    @cast_values.setter
    def cast_values(self, cast:bool):
        if cast:
            self._pairs = {
                k: _try_cast_to_numeric(_strip_quotes(v)) for (k,v) in self.pairs
            }
        else:
            warnings.warn(
                "Setting cast_values True->False has no effect on stored data.",
                category=ParseWarning, 
                stacklevel=2
            )
        self._cast_values = cast

    @property
    def pairs(self):
        """A dict containing key-value pairs extracted from the file.

        Numeric values will be parsed to int or float if possible. In these cases,
        precision specifiers will be stripped.

        Returns:
        --------
        dict[str : str | float | int]
        """
        return self._pairs

    @property
    def tables(self):
        """A list of data tables extracted from the file.

        These are stored as tuples of (keys, data) where data is a numpy array or list.

        Returns:
        --------
        list[(list[str], :math:`(N, N_{keys})` :class:`numpy.ndarray[str]`)]:
            A list of tuples corresponding with the keys and data of the system.
        """
        return self._tables


    @property
    def table_keys(self):
        """Column labels for tables extracted from the file.

        This is *only the labels* from the :meth:`~.tables` property.

        Returns:
        --------
        list[:math:`(N, N_{keys})` :class:`numpy.ndarray[str]`]:
            A list of arrays corresponding with the data of the system.
        """
        return [cols for (cols, data) in self.tables]


    @property
    def table_data(self):
        """Data for tables extracted from the file.

        This is *only the data* from the :meth:`~.tables` property.

        Returns:
        --------
        list[:math:`(N, N_{keys})` :class:`numpy.ndarray[str]`]:
            A list of arrays corresponding with the data of the system.
        """
        return [data for (cols, data) in self.tables]

    PATTERNS = {
        "key_value_numeric": r"^(_[\w\.]+)[ |\t]+(-?\d+\.?\d*)",
        "key_value_general": r"^(_[\w\.-]+)[ |\t]+([^#^\n]+)",
        "table_delimiter": r"([Ll][Oo][Oo][Pp]_)[ |\t]*([^\n]*)",
        "block_delimiter": r"([Dd][Aa][Tt][Aa]_)[ |\t]*([^\n]*)",
        "key_list": r"_[\w_\.]+",
        "space_delimited_data": r"(\'[^\']*\'|\"[^\"]*\"]|[^\'\" \t]*)[ | \t]*",
    }

    def __getitem__(self, key: str | list[str]):
        """Return an item from the dictionary of key-value pairs.

        Indexing with a string returns the value from the :meth:`~.pairs` dict. Indexing
        with an Iterable of strings returns a list of values, with `None` as a
        placeholder for keys that did not match any data.
        """
        if isinstance(key, Iterable) and not isinstance(key, str):
            return [self.pairs.get(k, None) for k in key]

        return self.pairs[key]
        # {key: data[key] for key in keys_to_include if key in data}

    def _parse(self):
        """Parse the cif file into python objects."""
        data_iter = peekable(self._data.split("\n"))

        for line in data_iter:
            if data_iter.peek(None) is None:
                break  # Exit without StopIteration

            # Combine nonsimple data values into a single, parseable line ==============
            while _line_is_continued(data_iter.peek()):
                line += _strip_comments(next(data_iter))
            line = _semicolon_to_string(line)

            # Skip processing if the line contains no data =============================
            if line == "" or _strip_comments(line) == "":
                continue

            # TODO: could support multi-block files in the future ======================
            block = re.match(self._cpat["block_delimiter"], line)
            if block is not None:
                continue

            # Extract key-value pairs and save to the internal state ===================
            pair = self._cpat["key_value_general"].match(line)
            if pair is not None:
                self._pairs.update(
                    {
                        pair.groups()[0]: _try_cast_to_numeric(
                            _strip_quotes(pair.groups()[1])
                        ) if self.cast_values else pair.groups()[1].strip()
                    }
                )

            # Build up tables by incrementing through the iterator =====================
            table = re.match(self._cpat["table_delimiter"], line)

            if table is not None:
                table_keys, table_data = [], []

                # First, extract table headers. Must be prefixed with underscore
                line_groups = table.groups()
                if line_groups[-1] != "":  # Extract table keys from the _loop line
                    fragment = _strip_comments(line_groups[-1].strip())
                    if fragment[:1] == "_":
                        keys = self._cpat["key_list"].findall(fragment)
                        table_keys.extend([] if keys is None else keys)
                    else:
                        continue

                while _is_key(data_iter.peek(None)):
                    line = _strip_comments(next(data_iter))
                    while _line_is_continued(data_iter.peek(None)):
                        line += _strip_comments(next(data_iter))
                    line = _semicolon_to_string(line)
                    table_keys.extend(self._cpat["key_list"].findall(line))

                while _is_data(data_iter.peek(None)):
                    line = _strip_comments(next(data_iter))
                    while _line_is_continued(data_iter.peek(None)):
                        line += _strip_comments(next(data_iter))
                    line = _semicolon_to_string(line)
                    parsed_line = self._cpat["space_delimited_data"].findall(line)
                    parsed_line = [_strip_quotes(m) for m in parsed_line if m != ""]
                    table_data.extend([parsed_line] if parsed_line else [])

                n_elements, n_cols = (
                    sum(len(row) for row in table_data),
                    len(table_keys),
                )

                if n_cols == 0:
                    continue  # Skip empty tables

                if n_elements % n_cols != 0:
                    warnings.warn(
                        f"Parsed data for table {len(self.tables)+1} cannot be resolved"
                        f" into a table of the expected size and will be returned as a "
                        f"list (got n={n_elements} items, expected c={n_cols} columns: "
                        f"n%c={n_elements % n_cols}).",
                        category=ParseWarning,
                        stacklevel=2,
                    )
                    # self._table_keys.append(table_keys)
                    # self._table_data.append(table_data)
                    self.tables.append((table_keys, table_data))
                    continue
                if not all(len(key) == len(table_keys[0]) for key in table_keys):
                    # warnings.warn(
                    #     f"Data for table {len(self.tables)+1} was parsed as a ragged arr",
                    #     category=ParseWarning,
                    #     stacklevel=1,
                    # )
                    table_data = np.array([*flatten(table_data)]).reshape(-1, n_cols)
                dt = _dtype_from_int(max(max(len(s) for s in l) for l in table_data))

                print("KEYS:", table_keys, f"N = {len(table_keys)}")
                print("DATA:", table_data)
                # print(max(max(len(s) for s in l) for l in table_data))
                print(dt)
                x = np.atleast_2d(table_data)
                x.dtype = [*zip(table_keys, [dt]*n_cols)]
                print(x)
                print(x.dtype)
                print(x["_publ_author_name"])
                self.tables.append((table_keys, np.atleast_2d(table_data)))

            if data_iter.peek(None) is None:
                break


if __name__ == "__main__":
    fn = "tests/sample_data/B-IncStrDb_Ccmm.cif"
    gen = _parsed_line_generator(fn, regexp=".*")

    cf = CifFile(fn=fn)
    # [print(pair) for pair in cf.pairs.items()]

    # print(cf.tables[0])
