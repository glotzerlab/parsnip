# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

import doctest

from sybil import Sybil
from sybil.parsers.doctest import DocTestParser

DOCTEST_OPTIONFLAGS = (
    doctest.NORMALIZE_WHITESPACE | doctest.IGNORE_EXCEPTION_DETAIL | doctest.ELLIPSIS
)


pytest_collect_file = Sybil(
    parsers=[
        # DocTestDirectiveParser(optionflags=DOCTEST_OPTIONFLAGS)
        DocTestParser(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS)
    ],
    pattern="*.rst",
    path="doc/source",
    # fixtures=['tempdir']
).pytest()
