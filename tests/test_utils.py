from pathlib import Path

import pytest

from parsnip._errors import ParseError, ParseWarning, _is_potentially_valid_path


def test_parse_error(capfd):
    with pytest.raises(ParseError) as error:
        raise ParseError("TEST_ERROR_RAISED")
    assert "TEST_ERROR_RAISED" in str(error.value)


def test_parse_warning():
    with pytest.raises(ParseWarning) as warning:
        raise ParseWarning("TEST_WARNING_RAISED")

    assert "TEST_WARNING_RAISED" in str(warning.value)


@pytest.mark.parametrize(
    ("path_str", "expected"),
    [
        (str(Path(__file__)), True),  # existing file
        (str(Path(__file__).parent / "conftest.py"), True),  # real file
        (str(Path(__file__).parent / "nonexistent.txt"), True),  # parent dir exists
        (str(Path(__file__).parent / "fake_file.cif"), True),  # .cif suffix
        (str(Path(__file__).parent / "asdf/noparent.txt"), False),  # no parent
        ("asdfasdfasd", False),
        ("asdfasdfasd.cif", True),
    ],
)
def test_is_potentially_valid_path(path_str, expected):
    assert _is_potentially_valid_path(path_str) is expected
