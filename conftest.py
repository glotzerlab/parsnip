# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.

import doctest

import pytest
from sybil import Sybil
from sybil.parsers.doctest import DocTestParser
from sybil.parsers.rest import DocTestDirectiveParser

DOCTEST_OPTIONFLAGS = (
    doctest.NORMALIZE_WHITESPACE | doctest.IGNORE_EXCEPTION_DETAIL | doctest.ELLIPSIS
)

@pytest.fixture(scope="module")
def changedir():
    import os
    cwd = os.getcwd()
    try:
        os.chdir("doc/source")
    finally:
        os.chdir(cwd)

pytest_collect_file = Sybil(
    parsers=[
         DocTestParser(optionflags=DOCTEST_OPTIONFLAGS),
         # DocTestDirectiveParser(optionflags=DOCTEST_OPTIONFLAGS)
        # DocTestParser(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS)
    ],
    pattern="*.rst",
    path="doc/source",
    fixtures=['changedir']
).pytest()
