# SPDX-License-Identifier: MIT
"""
Regression tests for Phase 3: state-mutating routes must reject requests
carrying an untrusted Origin header (CSRF defense), while GET requests and
requests with no Origin header at all (non-browser callers) stay unaffected.

Also guards against regressing the origin-regex anchoring bug found while
building this: the original CORS `allow_origin_regex` lacked a trailing '$'
on the localhost branch, so "http://localhost.evil.com" would have matched.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from fastapi.testclient import TestClient

from core.security import TRUSTED_ORIGIN_REGEX
from main import app


class TestRequireTrustedOriginOnRoutes(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)

    def test_post_with_no_origin_header_is_allowed(self):
        resp = self.client.post("/cache/flush")
        self.assertNotEqual(resp.status_code, 403)

    def test_post_with_trusted_localhost_origin_is_allowed(self):
        resp = self.client.post("/cache/flush", headers={"Origin": "http://localhost:1420"})
        self.assertNotEqual(resp.status_code, 403)

    def test_post_with_trusted_tauri_origin_is_allowed(self):
        resp = self.client.post("/cache/flush", headers={"Origin": "tauri://localhost"})
        self.assertNotEqual(resp.status_code, 403)

    def test_post_with_untrusted_origin_is_rejected(self):
        resp = self.client.post("/cache/flush", headers={"Origin": "https://evil.com"})
        self.assertEqual(resp.status_code, 403)

    def test_post_with_origin_suffix_attack_is_rejected(self):
        # Regression guard for the unanchored-regex bug: an attacker-owned
        # domain that merely starts with "localhost" must not pass.
        resp = self.client.post("/cache/flush", headers={"Origin": "http://localhost.evil.com"})
        self.assertEqual(resp.status_code, 403)

    def test_get_request_is_never_blocked_regardless_of_origin(self):
        resp = self.client.get("/jobs", headers={"Origin": "https://evil.com"})
        self.assertNotEqual(resp.status_code, 403)

    def test_delete_route_enforces_origin_check(self):
        resp = self.client.delete("/jobs/some-job-id", headers={"Origin": "https://evil.com"})
        self.assertEqual(resp.status_code, 403)

    def test_put_route_enforces_origin_check(self):
        resp = self.client.put(
            "/jobs/some-job-id/rename",
            json={"name": "x"},
            headers={"Origin": "https://evil.com"},
        )
        self.assertEqual(resp.status_code, 403)


class TestTrustedOriginRegexAnchoring(unittest.TestCase):
    """Unit-level guard on the shared regex itself, independent of routing."""

    def setUp(self):
        import re
        self.pattern = re.compile(TRUSTED_ORIGIN_REGEX)

    def test_legit_origins_match(self):
        for origin in [
            "http://localhost",
            "http://localhost:4200",
            "https://localhost:4200",
            "http://127.0.0.1",
            "http://127.0.0.1:8000",
            "tauri://localhost",
        ]:
            self.assertTrue(self.pattern.match(origin), f"expected {origin!r} to match")

    def test_suffix_and_subdomain_attacks_do_not_match(self):
        for origin in [
            "http://localhost.evil.com",
            "http://localhost-evil.com",
            "http://127.0.0.1.evil.com",
            "https://evil.com",
            "http://evil.com/localhost",
        ]:
            self.assertFalse(self.pattern.match(origin), f"expected {origin!r} to NOT match")


if __name__ == "__main__":
    unittest.main()
