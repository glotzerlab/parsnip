import numpy as np
import pytest

from parsnip._errors import ParseError, ParseWarning
from parsnip._utils import _deg2rad, _str2num


def test_parse_error():
    with pytest.raises(ParseError) as error:
        raise ParseError("TEST_ERROR_RAISED")

    assert "TEST_ERROR_RAISED" in str(error.value)


def test_parse_warning():
    with pytest.raises(ParseWarning) as warning:
        raise ParseWarning("TEST_WARNING_RAISED")

    assert "TEST_WARNING_RAISED" in str(warning.value)


def test_deg2rad(seed=43):
    rng = np.random.default_rng(seed)
    angles = rng.uniform(low=0, high=180, size=10_000)
    np.testing.assert_allclose(
        np.deg2rad(angles), [_deg2rad(val) for val in angles], atol=2e-15
    )


@pytest.mark.parametrize("string", ["3.1415926", "-12345", str(1e6), "0.00000003579"])
def test_str2num(string):
    converted_val = _str2num(string)
    if "." in string:
        assert isinstance(converted_val, float)
    else:
        assert isinstance(converted_val, int)
    assert np.isclose(float(string), converted_val)
