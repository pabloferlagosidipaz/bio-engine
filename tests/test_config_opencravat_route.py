# SPDX-License-Identifier: MIT
"""
Regression tests for Phase 4: /config/opencravat must validate that a
submitted oc_path resolves to a real executable before trusting it as
argv[0] of a later subprocess.run call (see services/opencravat.resolve_oc_executable).
"""
import os
import shutil
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from fastapi.testclient import TestClient

from core.config import settings
from main import app


class TestConfigureOpencravatValidation(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)
        self._original_oc_path = settings.oc_path

    def tearDown(self):
        settings.oc_path = self._original_oc_path

    def test_nonexistent_binary_is_rejected(self):
        resp = self.client.post(
            "/config/opencravat",
            json={"oc_path": "/definitely/not/a/real/binary/oc"},
        )
        self.assertEqual(resp.status_code, 400)
        # settings must not have been mutated by a rejected value
        self.assertEqual(settings.oc_path, self._original_oc_path)

    def test_resolvable_binary_is_accepted(self):
        # Use a real executable guaranteed to exist rather than "oc" itself,
        # so this test doesn't depend on OpenCRAVAT being installed.
        real_binary = shutil.which("python3") or shutil.which("python")
        self.assertIsNotNone(real_binary, "expected a real python executable on PATH for this test")

        resp = self.client.post("/config/opencravat", json={"oc_path": real_binary})
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(resp.json()["oc_path"], real_binary)
        self.assertEqual(settings.oc_path, real_binary)

    def test_empty_string_resets_to_default_oc_and_is_validated(self):
        # Empty string means "reset to default 'oc'" - if 'oc' isn't on PATH
        # in this environment, that must be surfaced as a 400, not silently
        # accepted as a broken default.
        oc_resolves = shutil.which("oc") is not None
        resp = self.client.post("/config/opencravat", json={"oc_path": ""})
        if oc_resolves:
            self.assertEqual(resp.status_code, 200)
            self.assertEqual(settings.oc_path, "oc")
        else:
            self.assertEqual(resp.status_code, 400)
            self.assertEqual(settings.oc_path, self._original_oc_path)

    def test_omitted_oc_path_leaves_setting_untouched(self):
        settings.oc_path = "some-preexisting-value"
        resp = self.client.post("/config/opencravat", json={})
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(settings.oc_path, "some-preexisting-value")


if __name__ == "__main__":
    unittest.main()
