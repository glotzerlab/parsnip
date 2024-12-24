import pytest

from parsnip.patterns import LineCleaner, _is_key, _is_data, _strip_comments, _strip_quotes, _semicolon_to_string
from parsnip._errors import ParseWarning

TEST_CASES = [
    None, "_key", "__key", "_key.loop_", "asdf", "loop_", "", " ", 
    "# comment", "_key#comment_ # 2", "loop_##com", "'my quote' # abc", 
    '"malformed \'\'#', ";oddness\"'\n;asdf", "_key.loop.inner_", 
    "normal_case", "multi.periods....", "__underscore__", 
    "_key_with_numbers123", "test#hash", "#standalone", 
    "'quote_in_single'", '"quote_in_double"', ' "mismatched_quotes\' ', 
    ";semicolon_in_text", ";;double_semicolon", "trailing_space ", 
    " leading_space", "_key.with#hash.loop", "__double#hash#inside__", 
    "single'; quote", "double;\"quote", "#comment;inside", 
    "_tricky'combination;#", ";'#another#combo;'", '"#edge_case"', 
    'loop;"#complex"case', "_'_weird_key'_", 
    "semi;;colon_and_hash#", "_odd.key_with#hash", 
    "__leading_double_underscore__", "middle;;semicolon", 
    "#just_a_comment", '"escaped \"quote"', 
    "'single_quote_with_hash#'", "_period_end.", "loop_.trailing_", 
    "escaped\\nnewline", "#escaped\\twith_tab", "only;semicolon", 
    "trailing_semicolon;", "leading_semicolon;", "_key;.semicolon", 
    "semicolon;hash#", "complex\"';hash#loop", "just_text", 
    "loop#weird\"text;", "nested'quotes\"here", "normal_case2", 
    "__underscored_case__", "escaped\\\"quotes#", ";semicolon#hash;", 
    "double#hash_inside##key", "__double..periods__", 
    "key#comment ; and_more", "_weird_;;#combo"
]


@pytest.mark.parametrize("line",TEST_CASES)
def test_is_key(line):
    if line is None or len(line)==0 or line[0] != "_":
        assert not _is_key(line)
        return
    assert _is_key(line)

@pytest.mark.parametrize("line",TEST_CASES)
def test_is_data(line):
    if line is not None and len(line) == 0:
        assert _is_data(line)
    elif line is None or line[0] == "_" or line[:5] == "loop_":
        assert not _is_data(line)
        return
    assert _is_data(line)

@pytest.mark.parametrize("line",TEST_CASES)
def test_strip_comments(line):
    if line is None:
        return
    elif all(c==" " for c in line):
        assert _strip_comments(line) == ""
        return
    elif "#" not in line and not all(c==" " for c in line):
        assert _strip_comments(line) == line.strip()
        return

    stripped = _strip_comments(line)
    assert "#" not in stripped
    assert len(stripped) < len(line)

@pytest.mark.parametrize("line",TEST_CASES)
def test_strip_quotes(line):
    if line is None:
        return
    elif "'" not in line and '"' not in line:
        assert _strip_quotes(line) == line
        return

    stripped = _strip_quotes(line)
    assert "'" not in stripped and '"' not in stripped
    assert len(stripped) < len(line)

@pytest.mark.parametrize(
    "line",
    # [None, "_key", "__key", "_key.loop_", "asdf", "loop_", "", " ", "# comment", "_key#comment_ # 2", "loop_##com", "'my quote' # abc", '"malformed \'\'#', ";oddness\"'\n;asdf"]
    TEST_CASES
)
def test_semicolon_to_string(line):
    if line is None:
        return
    elif "'" in line and '"' in line:
        with pytest.warns(ParseWarning, match="String contains single and double quotes"):
            fixed = _semicolon_to_string(line)
            assert (fixed == line) if ";" not in line else (";" not in fixed)
            return
    elif ";" not in line:
        assert _semicolon_to_string(line) == line
        return
    fixed = _semicolon_to_string(line)
    assert ";" not in fixed
    assert "'" in fixed if "'" not in line else '"' in fixed

