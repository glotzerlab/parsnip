"""Functions and classes to process string data."""
import re

# Compile in common patterns for cif parsing. These are reused throughout the package.
_multiple_whitespace_pattern = re.compile(r"\s+")
_comma_prune_spaces = re.compile(r",\s+")


def compile_pattern_from_strings(filter_patterns: tuple[str]):
    """Return a regex pattern that matches any of the characters in the filter.

    Args:
        filter_patterns (list[str]): Description

    Returns:
        re.Pattern: Pattern matching any of the input characters.
    """
    return re.compile("|".join(filter_patterns))


class LineCleaner:
    """Simple object to apply a series of regex patterns to a string.

    To intialize a line cleaner, pass in a tuple of strings of the form
    ``(pattern, replacement)``. Patterns are compiled on initialization to accelerate
    future processing.

    Args:
        patterns (tuple[tuple[str,str]]): Tuple of tuples of strings.
            The first item in each tuple is the pattern to match, and the second item is
            what that pattern will be replaced with.
    """

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

    def __call__(self, line: str):
        """Apply patterns defined on initialization of the object to the string.

        ``re.sub(pattern,line)`` is run for each pattern (in order) in self.patterns,
        which is defined on initialization.

        Args:
            line (str): String to apply patterns to.

        Returns:
            str: The substituted lines.
        """
        for pattern, replacement in zip(self.patterns, self.replacements):
            line = pattern.sub(replacement, line)
        return line
