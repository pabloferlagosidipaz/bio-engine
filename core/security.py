# SPDX-License-Identifier: MIT
"""
CSRF defense for state-mutating routes.

CORS (see main.py) only controls whether a cross-origin script can *read* a
response - it does not stop a browser from *sending* the request in the
first place. A multipart/form-data POST is a CORS "simple request" (no
preflight), so any webpage the user happens to have open can silently
trigger routes like /upload against this local server. Browsers attach an
Origin header on every cross-origin POST/PUT/DELETE regardless of CORS, so
we validate it server-side here rather than relying on CORS alone.
"""
import re

from fastapi import Request

from core.exceptions import BioEngineError

# Single source of truth for "who is allowed to talk to this sidecar" - also
# used by main.py for the CORS allow_origin_regex. Each alternative is fully
# anchored (^...$); an earlier version of this pattern lacked the trailing
# '$' on the localhost branch, which meant an origin like
# "http://localhost.evil.com" (an attacker-registered domain) would also
# match, since re.match only anchors the start of the string.
TRUSTED_ORIGIN_REGEX = (
    r"^https?://localhost(:\d+)?$"
    r"|^https?://127\.0\.0\.1(:\d+)?$"
    r"|^tauri://localhost$"
)

_trusted_origin_pattern = re.compile(TRUSTED_ORIGIN_REGEX)
_UNSAFE_METHODS = {"POST", "PUT", "DELETE", "PATCH"}


def require_trusted_origin(request: Request) -> None:
    """
    FastAPI dependency: rejects state-mutating requests whose Origin header
    (when present) isn't in the trusted set defined by TRUSTED_ORIGIN_REGEX.

    Requests with no Origin header at all (curl, most native Tauri IPC
    calls, test clients) are allowed through - this specifically closes the
    "malicious webpage open in the user's browser" CSRF vector, not general
    local-process access, which remains an accepted local-trust assumption
    (see code review report).
    """
    if request.method not in _UNSAFE_METHODS:
        return

    origin = request.headers.get("origin")
    if origin is None:
        return

    if not _trusted_origin_pattern.match(origin):
        raise BioEngineError(
            f"Requests from origin '{origin}' are not permitted",
            context="csrf_protection",
            status_code=403,
        )
