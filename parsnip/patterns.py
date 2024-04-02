"""TODO."""
import re

# Compile in common patterns for cif parsing. These are reused throughout the package.
_multiple_whitespace_pattern = re.compile(r"\s+")
_comma_prune_spaces = re.compile(r",\s+")


def compile_pattern_from_strings(strings: tuple[str]):
    """Return a regex pattern that matches any of the characters in the filter.

    Args:
        strings (list[str]): Description

    Returns:
        re.Pattern: Pattern matching any of the input characters.
    """
    return re.compile("|".join(strings))


class LineCleaner:
    """TODO."""

    def __init__(self, patterns: tuple[tuple[str, str]]):
        self.patterns, self.replacements = [], []

        # If we only have a single tuple
        if isinstance(patterns[0], str):
            pattern, replacement = patterns
            self.patterns.append(re.compile(pattern))

            self.replacements.append(replacement)
        else:
            for pattern, replacement in patterns:
                self.patterns.append(re.compile(pattern))

                self.replacements.append(replacement)

    def __call__(self, line):
        """A."""
        for pattern, replacement in zip(self.patterns, self.replacements):
            # print(pattern,replacement)
            line = pattern.sub(replacement, line)
        return line
