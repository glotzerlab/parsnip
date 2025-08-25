import numpy as np
import pytest
from conftest import cif_files_mark

from parsnip._errors import ParseWarning
from parsnip.patterns import (
    _box_from_lengths_and_angles,
    _dtype_from_int,
    _is_data,
    _is_key,
    _semicolon_to_string,
    _strip_comments,
    _strip_quotes,
    _try_cast_to_numeric,
    _write_debug_output,
)

BOX_ATOL = 8e-7

TEST_CASES = [
    None,
    "_key",
    "__key",
    "_key.loop_",
    "asdf",
    "loop_",
    "",
    " ",
    "# comment",
    "_key#comment_ # 2",
    "loop_##com",
    "'my quote' # abc",
    "\"malformed ''#",
    ";oddness\"'\n;asdf",
    "_key.loop.inner_",
    "normal_case",
    "multi.periods....",
    "__underscore__",
    "_key_with_numbers123",
    "test#hash",
    "#standalone",
    "'quote_in_single'",
    '"quote_in_double"',
    " \"mismatched_quotes' ",
    ";semicolon_in_text",
    ";;double_semicolon",
    "trailing_space ",
    " leading_space",
    "_key.with#hash.loop",
    "__double#hash#inside__",
    "single'; quote",
    'double;"quote',
    "#comment;inside",
    "_tricky'combination;#",
    ";'#another#combo;'",
    '"#edge_case"',
    'loop;"#complex"case',
    "_'_weird_key'_",
    "semi;;colon_and_hash#",
    "_odd.key_with#hash",
    "__leading_double_underscore__",
    "middle;;semicolon",
    "#just_a_comment",
    '"escaped "quote"',
    "'single_quote_with_hash#'",
    "_period_end.",
    "loop_.trailing_",
    "escaped\\nnewline",
    "#escaped\\twith_tab",
    "only;semicolon",
    "trailing_semicolon;",
    "leading_semicolon;",
    "_key;.semicolon",
    "semicolon;hash#",
    "complex\"';hash#loop",
    "just_text",
    'loop#weird"text;',
    "nested'quotes\"here",
    "normal_case2",
    "__underscored_case__",
    'escaped\\"quotes#',
    ";semicolon#hash;",
    "double#hash_inside##key",
    "__double..periods__",
    "key#comment ; and_more",
    "_weird_;;#combo",
]


@pytest.mark.parametrize("line", TEST_CASES)
def test_is_key(line):
    if line is None or len(line) == 0 or line[0] != "_":
        assert not _is_key(line)
        return
    assert _is_key(line)


@pytest.mark.parametrize("line", TEST_CASES)
def test_is_data(line):
    if line is not None and len(line) == 0:
        assert _is_data(line)
    elif line is None or line[0] == "_" or line[:5] == "loop_":
        assert not _is_data(line)
        return
    assert _is_data(line)


@pytest.mark.parametrize("line", TEST_CASES)
def test_strip_comments(line):
    if line is None:
        return
    if all(c == " " for c in line):
        assert _strip_comments(line) == ""
        return
    if "#" not in line and not all(c == " " for c in line):
        assert _strip_comments(line) == line.strip()
        return

    stripped = _strip_comments(line)
    assert "#" not in stripped
    assert len(stripped) < len(line)


@pytest.mark.parametrize("line", TEST_CASES)
def test_strip_quotes(line):
    if line is None:
        return
    if "'" not in line and '"' not in line:
        assert _strip_quotes(line) == line
        return

    stripped = _strip_quotes(line)
    assert "'" not in stripped
    assert '"' not in stripped
    assert len(stripped) < len(line)


@pytest.mark.parametrize("line", TEST_CASES)
def test_semicolon_to_string(line):
    if line is None:
        return
    if "'" in line and '"' in line:
        with pytest.warns(
            ParseWarning, match="String contains single and double quotes"
        ):
            fixed = _semicolon_to_string(line)
        assert (fixed == line) if ";" not in line else (";" not in fixed)
        return
    if ";" not in line:
        assert _semicolon_to_string(line) == line
        return
    fixed = _semicolon_to_string(line)
    assert ";" not in fixed
    assert "'" in fixed if "'" not in line else '"' in fixed


@pytest.mark.parametrize("nchars", range(1, 15))
def test_dtype_from_int(nchars):
    assert _dtype_from_int(nchars) == "<U" + str(nchars)


def has_duplicates(pos):
    pos = np.array(pos)
    distances = np.linalg.norm(pos[:, None] - pos, axis=-1)
    np.fill_diagonal(distances, np.inf)  # Ignore self-comparison
    unique_indices = np.arange(len(pos))[distances[:, 0] != 0]
    unique_counts = np.sum(distances == 0, axis=1)[unique_indices] + 1

    return unique_indices, unique_counts


@pytest.mark.parametrize(
    "pos",
    [
        [[1, 2, 3], [3, 4, 5], [6, 7, 8], [9, 10, 11]],
        [[0, 1, 2], [0, 1, 2], [3, 4, 5], [6, 7, 8]],
        [[0, 1, 2], [0, 1, 2], [3, 4, 5], [3, 4, 5]],
        [[0, 1, 2], [0, 1, 2], [0, 1, 2], [6, 7, 8]],
        [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]],
    ],
)
def test_write_debug_output(pos, capfd):
    computed_indices, computed_counts = has_duplicates(pos)
    _write_debug_output(computed_indices, computed_counts, pos)

    duplicates_exist = any(count > 1 for count in computed_counts)
    if duplicates_exist:
        assert "(duplicate point, number of occurrences)" in capfd.readouterr().out
    else:
        assert "... all points are unique (within tolerance)." in capfd.readouterr().out


@pytest.mark.parametrize(
    "s", ["1.234", "abcd", "1999", "33(45)", "01.2", "8.9(1)", "9.87a"]
)
def test_try_cast_to_numeric(s):
    result = _try_cast_to_numeric(s)

    if "a" in s:
        assert isinstance(result, str)
    elif "." in s:
        assert isinstance(result, float)
    else:
        assert isinstance(result, int)


@cif_files_mark
def test_repr(cif_data):
    import re

    repr = re.sub(r"[a-z\s]*", "", cif_data.file.__repr__().split(" : ")[1]).split(",")
    n_pairs, n_tabs = [int(i) for i in repr]
    assert n_pairs == len(cif_data.file.pairs)
    assert n_tabs == len(cif_data.file.loops)


@cif_files_mark
def test_box(cif_data):
    freud = pytest.importorskip("freud")

    cif_box = cif_data.file.read_cell_params(degrees=False)

    freud_box = freud.Box.from_box_lengths_and_angles(*cif_box)
    freud_box_2 = freud.Box(*cif_data.file.box)
    parsnip_box = _box_from_lengths_and_angles(*cif_box)

    np.testing.assert_allclose(parsnip_box[:3], freud_box.L, atol=BOX_ATOL)
    np.testing.assert_allclose(
        parsnip_box[3:], [freud_box.xy, freud_box.xz, freud_box.yz], atol=BOX_ATOL
    )
    if "PDB" not in cif_data.filename:
        np.testing.assert_allclose(
            [*freud_box.L, freud_box.xy, freud_box.xz, freud_box.yz],
            [*freud_box_2.L, freud_box_2.xy, freud_box_2.xz, freud_box_2.yz],
            atol=BOX_ATOL,
        )
