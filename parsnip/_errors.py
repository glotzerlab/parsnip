# Copyright (c) 2025, The Regents of the University of Michigan
# This file is from the parsnip project, released under the BSD 3-Clause License.

import warnings
from pathlib import Path


def _is_potentially_valid_path(file: str) -> bool:
    """Check whether a file string could possibly be intended as a path.

    This method returns true if the provided string is a valid path, whether the suffix
    ".cif" is contained in the path, if the path links to a file, or if the path's
    parent is a directory.
    """
    try:
        path = Path(file)
        return (
            ".cif" in path.suffixes  # Probably intended to parse as file
            or path.exists()  # If it is a file, we definitely want to parse that
            # Possibly a typo, but we want to check that path regardless.
            or (path.parent.is_dir() and path.parent != Path("."))
        )
    except OSError:
        return False


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


def _warn_or_err(msg, strict):
    if strict:
        raise ValueError(msg)
    warnings.warn(
        msg.replace("\n", ""),
        category=ParseWarning,
        stacklevel=2,
    )
