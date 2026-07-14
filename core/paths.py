# SPDX-License-Identifier: MIT
"""
Filesystem path safety helpers.

Any path fragment that originates from a client request (an upload filename,
a query-string path, ...) must never be trusted directly in a filesystem
operation - see CWE-22 (Path Traversal). Route handlers that accept such
input MUST go through one of the helpers below rather than calling
os.path.join / open() directly on client-supplied strings.
"""
import os

from core.exceptions import BioEngineError


class PathTraversalError(BioEngineError):
    """Raised when a client-supplied path/filename escapes its allowed directory."""
    def __init__(self, message: str = "Invalid path"):
        super().__init__(message, context="path_validation", status_code=400)


def safe_filename_join(base_dir: str | os.PathLike, filename: str) -> str:
    """
    Joins a client-supplied *filename* (not a path) onto a trusted base
    directory, guaranteeing the result stays inside base_dir.

    Only the basename of `filename` is used, so any directory components
    (including `..` or a leading `/`) are stripped before joining. Use this
    for inputs that are conceptually "just a filename", e.g. an uploaded
    file's original name.
    """
    if not filename:
        raise PathTraversalError("Filename must not be empty")

    name = os.path.basename(filename)
    if not name or name in (".", ".."):
        raise PathTraversalError(f"Invalid filename: {filename!r}")

    base_real = os.path.realpath(base_dir)
    candidate = os.path.realpath(os.path.join(base_real, name))

    if os.path.commonpath([base_real, candidate]) != base_real:
        raise PathTraversalError(f"Invalid filename: {filename!r}")

    return candidate


def resolve_within(path: str, allowed_dirs: list[str | os.PathLike]) -> str:
    """
    Resolves a client-supplied *path* and guarantees it falls inside one of
    `allowed_dirs`. Use this for inputs that are conceptually a full path
    (e.g. a path previously handed back by /upload), where restricting to a
    bare filename isn't appropriate.

    Returns the resolved (realpath) string. Raises PathTraversalError if the
    path doesn't resolve inside any allowed directory.
    """
    if not path:
        raise PathTraversalError("Path must not be empty")

    candidate = os.path.realpath(path)

    for allowed in allowed_dirs:
        allowed_real = os.path.realpath(allowed)
        if os.path.commonpath([allowed_real, candidate]) == allowed_real:
            return candidate

    raise PathTraversalError(f"Path is outside of allowed directories: {path!r}")
