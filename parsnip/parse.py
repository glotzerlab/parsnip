"""CIF parsing tools."""

import warnings

import numpy as np

from ._utils import ParseError, ParseWarning
from .patterns import LineCleaner


def _remove_comments_from_line(line):
    return line.split("#")[0].strip()


def read_table(
    filename: str,
    keys: str,
    filter_line: tuple[tuple[str, str]] = ((r",\s+", ",")),
    keep_original_key_order=False,
) -> np.ndarray[str]:
    r"""Extract data from a CIF file loop_ table.

    CIF files store tabular data as whitespace-delimited blocks that start with `loop_`.
    Keys are kept at the top of the table, and the vertical position of keys corresponds
    to the horizontal position of the column storing the data for that key. The end of
    the table is not necessarily marked: instead, the script detects when the table
    format is exited.

    For example:

    ```
    loop_
    _space_group_symop_id
    _space_group_symop_operation_xyz
    1 x,y,z
    2 -x,y,-z+1/2
    3 -x,-y,-z
    4 x,-y,z+1/2
    5 x+1/2,y+1/2,z
    6 -x+1/2,y+1/2,-z+1/2
    7 -x+1/2,-y+1/2,-z
    8 x+1/2,-y+1/2,z+1/2

    ```

    Only data columns corresponding to a key in the input keys list will be returned.

    Note that this function will ONLY return data from a single table. If keys are
    provided that correspond to data from multiple tables, only the first table will
    be read.

    The ``filter_line`` argument allows for dynamic input creation of regex filters to
    apply to each line that contains data to be saved. The default value is
    ``((",\s+",","))``, which helps differentiate between individual data fragments
    seperated by commas and whitespace characters, and other sections of the line that
    are also whitespace separated. Adding another tuple to remove single quotes can
    also be helpful: try ``((",\s+",","),(",",""))`` to achieve this. To disable the
    feature entirely, pass in a tuple of empty strings: ``("","")``. Note that doing so
    will cause errors if the table contains non-delimiting whitespaces.

    Args:
        filename (str): The name of the .cif file to be parsed.
        keys (tuple[str]): The names of the keys to be parsed.
        filter_line (tuple[tuple[str]], optional):
            A tuple of strings that are compiled to a regex filter and applied to each
            data line. (Default value: ((r",\s+",",")) )
        keep_original_key_order (bool, optional):
            When True, preserve the order of keys in the table from the cif file.
            When False, return columns of data in order of the input ``keys`` arg.
            (Default value: False)

    Returns:
        np.ndarray[str]: A numpy array of the data as strings.
    """
    with open(filename) as f:
        tables = f.read().split("loop_")

    line_cleaner = LineCleaner(filter_line)
    nontable_line_prefixes = ("_", "#")

    for table in tables:
        lines = table.strip().split("\n")
        in_header = True
        data_column_indices, data, column_order = [], [], []

        for line_number, line in enumerate(lines):
            # Check for invalid blank lines in the table header
            if in_header and data_column_indices and line == "":
                raise ParseError(
                    "Whitespace may not be used in between keys in the table header. "
                    "See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#general"
                    ", section 7 for more details."
                )

            # We will get errors if there is a comment after the loop_ block that
            # contains our data. This is questionably legal, but very uncommon

            line = _remove_comments_from_line(line)

            # Save current key position if it is one of the keys we want.
            if in_header and (line in keys):
                data_column_indices.append(line_number)
                if not keep_original_key_order:
                    column_order.append(keys.index(line))
                continue

            # If we exit the header and enter the table body
            if data_column_indices and (line[:1] not in nontable_line_prefixes):
                in_header = False  # Exit the header and start writing data
                clean_line = line_cleaner(line)
                split_line = clean_line.split()

                # Only add data if the line has at least as many columns as required.
                n_cols_found, n_cols_expected = (
                    len(split_line),
                    len(data_column_indices),
                )
                if n_cols_found >= n_cols_expected:
                    data.append(split_line)
                elif split_line != [] and n_cols_found < n_cols_expected:
                    warnings.warn(
                        ParseWarning(
                            f"Data line is a fragment and will be skipped: (expected "
                            f"line with {n_cols_expected} values, got {split_line})."
                        ),
                        stacklevel=2,
                    )
                continue
            elif (not in_header) and (line[:1] == "_"):
                break
        if data_column_indices:
            break

    if not keep_original_key_order:
        # Reorder the column indices to match the order of the input keys
        data_column_indices = np.array(data_column_indices)[np.argsort(column_order)]

    return np.atleast_2d(data)[:, data_column_indices]
