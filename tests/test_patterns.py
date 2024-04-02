import pytest

from parsnip.patterns import LineCleaner


@pytest.mark.parametrize("string", ["x,y,      z", "'x,  y, z'"])
def test_linecleaner(string):
    pattern = ([r"\s+", " "], (r",\s+", ","))
    print(pattern)
    cleaner = LineCleaner(pattern)

    if string == "x,y,z" and pattern == ((r"\s+", ",")):
        return

    assert cleaner(string) == "x,y,z"
