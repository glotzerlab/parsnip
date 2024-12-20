# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

import re
import warnings

from more_itertools import peekable

from parsnip.parse import _parsed_line_generator

NONTABLE_LINE_PREFIXES = ("_", "#")


class CifFile:
    """Parsed CIF file."""

    def __init__(self, fn: str):
        self._fn = fn
        self._pairs = {}

        self._cpat = {k: re.compile(pattern) for (k, pattern) in self.PATTERNS.items()}

        with open(fn) as file:
            self._data = file.read()
        self._parse()

    @property
    def pairs(self):
        return self._pairs

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
            return data_iter.peek().strip()[:1] == ";"

        data_iter = peekable(self._data.split("\n"))

        for line in data_iter:
            if data_iter.peek(None) is None:
                break  # Exit without StopIteration

            # Combine nonsimple data values into a single, parseable line
            while line_is_continued(line):
                print("CONTINUING")
                line += strip_comments(next(data_iter))
            line = semicolon_to_string(line)
            if line == "" or strip_comments(line) == "":
                continue

            table = re.match(self._cpat["table_delimiter"], line)

            if table is not None:
                print("ENTERING TABLE")
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

                # Extract table headers from subsequent lines
                def is_key(line):
                    stripped = line.strip()
                    return stripped != "" and stripped[0] == "_"

                def is_data(line):
                    stripped = line.strip()
                    return stripped[:1] != "_"

                line = strip_comments(next(data_iter))

                while is_key(line):
                    table_keys.extend(self._cpat["key_list"].findall(line))
                    line = strip_comments(next(data_iter))

                while is_data(line):
                    parsed_line = self._cpat["space_delimited_data"].findall(line)
                    table_data.extend([parsed_line] if parsed_line else [])
                    line = strip_comments(next(data_iter))

                print("KEYS:", table_keys)
                print("VALS:", table_data)

            # TODO: could support multi-block files in the future
            block = re.match(self._cpat["block_delimiter"], line)
            if block is not None:
                continue

            print(line)
            print(self._cpat["key_value_general"].match(line).groups())
            print()
            # self.pairs.update({})

            if next(data_iter, None) is None:
                break


fn = "tests/sample_data/B-IncStrDb_Ccmm.cif"
gen = _parsed_line_generator(fn, regexp=".*")

cf = CifFile(fn=fn)
# print(cf._data)


# print([*gen])
# print()
# print([*gen])
