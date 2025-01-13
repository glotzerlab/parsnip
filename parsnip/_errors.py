# Copyright (c) 2024, Glotzer Group
# This file is from the parsnip project, released under the BSD 3-Clause License.


class ParseWarning(Warning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class ParseError(RuntimeError):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)
