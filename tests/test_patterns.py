import pytest

from parsnip.patterns import LineCleaner


@pytest.mark.parametrize("string", ["x,y,      z", "'x,  y, z'", "'x',  'y', 'z'"])
def test_linecleaner_xyz(string):
    patterns = ([r",\s+", ","], ("'", ""))
    cleaner = LineCleaner(patterns)

    assert cleaner(string) == "x,y,z"


@pytest.mark.parametrize(
    "string",
    [
        "As4 As+3 0.3172 0.9426 1.017 1 0.0",
        "O(1)  1   0.179(1)         0  0.791(4)  ?  Uiso",
        "Br02 Br 0.0000 0.5000 0.5000 0.0983(7) Uani 1 16 d S T P . .",
        "C1 C 0.0000 -0.120(5) 0.0000 0.030(12) Uiso 0.0417 8 d DS . P . 1",
    ],
)
def test_linecleaner_floatstrip(string):
    patterns = [r"\(\d*\)", ""]
    cleaner = LineCleaner(patterns)

    assert cleaner(string).split() == [
        substr.split("(")[0] for substr in string.split()
    ]
