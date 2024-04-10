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
