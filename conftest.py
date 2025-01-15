import doctest

from sybil import Sybil
from sybil.parsers.codeblock import PythonCodeBlockParser, CodeBlockParser
from sybil.parsers.doctest import DocTestParser
from sybil.parsers.rest import DocTestDirectiveParser

DOCTEST_OPTIONFLAGS = (
    doctest.NORMALIZE_WHITESPACE | doctest.IGNORE_EXCEPTION_DETAIL | doctest.ELLIPSIS
)

pytest_collect_file = Sybil(
    parsers=[
        DocTestDirectiveParser(optionflags=DOCTEST_OPTIONFLAGS)
    ],
    pattern='*.rst',
    # fixtures=['tempdir']
).pytest()
