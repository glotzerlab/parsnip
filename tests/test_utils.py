import pytest

from parsnip._errors import ParseError, ParseWarning


def test_parse_error(capfd):
    with pytest.raises(ParseError) as error:
        raise ParseError("TEST_ERROR_RAISED")
    assert "TEST_ERROR_RAISED" in str(error.value)


def test_parse_warning():
    with pytest.raises(ParseWarning) as warning:
        raise ParseWarning("TEST_WARNING_RAISED")

    assert "TEST_WARNING_RAISED" in str(warning.value)
