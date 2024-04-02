"""CIF parsing tools."""


import numpy as np

from ._utils import ParseWarning
from .patterns import LineCleaner


def read_table(
    filename: str,
    keys: str,
    filter_line: tuple[tuple[str, str]] = ((r",\s+", ",")),
    validate: bool = True,
) -> np.ndarray[str]:
    r"""Extract data from a CIF file loop_ table.

    CIF files store tabular data as whitespace-delimited blocks surrounded by `loop_`.
    Keys are kept at the top of the table, and the vertical position of keys corresponds
    to the horizontal position of the column storing the data for that key.
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

    loop_
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
    feature entirely, pass in a tuple of empty strings: ``("","")``

    Args:
        filename (str): The name of the .cif file to be parsed.
        keys (tuple[str]): The names of the keys to be parsed.
        filter_line (tuple[tuple[str]], optional):
            A tuple of strings that are compiled to a regex filter and applied to each
            data line.
            (Default value: ((r",\s+",",")) )

        validate (bool, optional):
            Whether to check the output is correct. (Default value: True)

    Returns:
        np.ndarray[str]: A numpy array of the data as strings.
    """
    with open(filename) as f:
        blocks = f.read().split("loop_")

    line_cleaner = LineCleaner(filter_line)

    for table in blocks:
        lines = table.strip().split("\n")
        in_header = True
        data_column_indices, data = [], []

        for line_number, line in enumerate(lines):
            if in_header and (line in keys):
                data_column_indices.append(line_number)
                continue

            # Take a slice to avoid indexing a 0 length string
            if data_column_indices and line[:1] != "_":
                in_header = False  # Exit the header and start writing data

                clean_line = line_cleaner(line)
                split_line = clean_line.split()

                # Check that we are not appending a fragment line
                if len(split_line) == len(data_column_indices):
                    data.append(split_line)
                # Warn if we get a non-empty line with the wrong number of values
                elif split_line != []:
                    raise ParseWarning(
                        f"Data line is a fragment: (expected line with "
                        f"{len(data_column_indices)} values, got {split_line})"
                    )

            elif data_column_indices and line[:1] == "_":
                break  # If we are back to standard data, exit the loop

        if data_column_indices:  # Exit the loop if we've found our data.
            n_keys, n_cols = len(keys), len(data_column_indices)
            if n_keys != n_cols:
                raise ParseWarning(
                    "Data has incorrect number of columns: "
                    f"(expected {n_keys}, got {n_cols})"
                )
            break

    return np.atleast_2d(data)[:, data_column_indices]
