# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

import re
import warnings

from more_itertools import peekable

from parsnip.parse import _parsed_line_generator

NONTABLE_LINE_PREFIXES = ("_", "#")

def is_key(line):
    return line is not None and line.strip()[:1] == "_"

def is_data(line):
    return line is not None and line.strip()[:1] != "_"

class CifFile:
    """Parsed CIF file."""

    def __init__(self, fn: str):
        self._fn = fn
        self._pairs = {}
        self._tables = []

        self._cpat = {k: re.compile(pattern) for (k, pattern) in self.PATTERNS.items()}

        with open(fn) as file:
            self._data = file.read()
        self._parse()

    @property
    def pairs(self):
        return self._pairs

    @property
    def tables(self):
        return self._tables

    PATTERNS = {
        "key_value_numeric": r"^(_[\w\.]+)[ |\t]+(-?\d+\.?\d*)",
        "key_value_general": r"^(_[\w\.-]+)[ |\t]+([^#^\n]+)",
        "table_delimiter": r"([Ll][Oo][Oo][Pp]_)[ |\t]*([^\n]*)",
        "block_delimiter": r"([Dd][Aa][Tt][Aa]_)[ |\t]*([^\n]*)",
        "key_list": r"_[\w_]+",
        "space_delimited_data": r"'[^']*'|\"[^\"]*\"|\S+",
    }

    def _parse(self):
        def strip_comments(s: str):
            return s.split("#")[0].strip()

        def semicolon_to_string(line: str):
            if "'" in line and '"' in line:
                warnings.warn(
                    (
                        "String contains single and double quotes - "
                        "line may be parsed incorrectly"
                    ),
                    stacklevel=2,
                )
            # WARNING: because we split our string, we strip "\n" implicitly
            return line.replace(";", "'" if "'" not in line else '"')

        def line_is_continued(line: str):
            return line.strip()[:1] == ";"

        data_iter = peekable(self._data.split("\n"))

        for line in data_iter:
            if data_iter.peek(None) is None:
                break  # Exit without StopIteration

            # Combine nonsimple data values into a single, parseable line ==============
            while line_is_continued(data_iter.peek()):
                line += strip_comments(next(data_iter))
            line = semicolon_to_string(line)

            # Skip processing if the line contains no data =============================
            if line == "" or strip_comments(line) == "":
                continue

            # TODO: could support multi-block files in the future ======================
            block = re.match(self._cpat["block_delimiter"], line)
            if block is not None: continue

            # Extract key-value pairs and save to the internal state ===================
            pair = self._cpat["key_value_general"].match(line)
            if pair is not None:
                self._pairs.update({pair.groups()[0]: pair.groups()[1]})


            # Build up tables by incrementing through the iterator =====================
            table = re.match(self._cpat["table_delimiter"], line)

            if table is not None:
                table_keys, table_data = [], []

                # First, extract table headers. Must be prefixed with underscore
                line_groups = table.groups()
                if line_groups[-1] != "":  # Extract table keys from the _loop line
                    fragment = strip_comments(line_groups[-1].strip())
                    if fragment[:1] == "_":
                        keys = self._cpat["key_list"].findall(fragment)
                        table_keys.extend([] if keys is None else keys)
                    else:
                        continue

                while is_key(data_iter.peek(None)):
                    line = strip_comments(next(data_iter))
                    table_keys.extend(self._cpat["key_list"].findall(line))

                while is_data(data_iter.peek(None)):
                    line = strip_comments(next(data_iter))
                    parsed_line = self._cpat["space_delimited_data"].findall(line)
                    table_data.extend([parsed_line] if parsed_line else [])

                print("KEYS:", table_keys)
                print("VALS:", table_data)


            if data_iter.peek(None) is None:
                break
        [print(pair) for pair in self._pairs.items()]

fn = "tests/sample_data/B-IncStrDb_Ccmm.cif"
gen = _parsed_line_generator(fn, regexp=".*")

cf = CifFile(fn=fn)
# print(cf._data)


# print([*gen])
# print()
# print([*gen])
