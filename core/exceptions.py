class BioEngineError(Exception):
    """Base exception for the bio-engine."""
    def __init__(self, message: str, context: str = "processing"):
        super().__init__(message)
        self.message = message
        self.context = context

class ReferenceError(BioEngineError):
    """Raised when there is an issue with the reference file or NCBI ID."""
    pass

class AlignmentError(BioEngineError):
    """Raised when alignment fails."""
    pass

class FileNotFoundError(BioEngineError):
    """Raised when a required file is missing."""
    pass
