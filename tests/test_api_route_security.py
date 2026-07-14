# SPDX-License-Identifier: MIT
"""
Security regression tests for path-handling API routes.

Covers:
  - Phase 1: /upload must not allow a crafted filename to write outside
    settings.uploads_dir (CWE-22 / arbitrary file write).
  - Phase 2: /preview-read must not allow an arbitrary `path` query param
    to read a file outside settings.uploads_dir / settings.jobs_dir
    (CWE-22 / arbitrary file read).

These tests exercise the real FastAPI app through TestClient so they catch
regressions in the route wiring itself, not just the helper functions in
core/paths.py (see test_paths.py for helper-level unit tests).
"""
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from fastapi.testclient import TestClient

from core.config import settings
from main import app

_FAKE_TRACE = {
    "traceA": [], "traceC": [], "traceG": [], "traceT": [],
    "peakLocations": [], "basecalls": {},
}


class TestUploadPathTraversal(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)
        self.tmp = tempfile.TemporaryDirectory()
        self.uploads_dir = os.path.join(self.tmp.name, "uploads")
        os.makedirs(self.uploads_dir, exist_ok=True)
        self._patcher = patch.object(settings, "uploads_dir", self.uploads_dir)
        self._patcher.start()

    def tearDown(self):
        self._patcher.stop()
        self.tmp.cleanup()

    def test_normal_filename_uploads_into_uploads_dir(self):
        resp = self.client.post(
            "/upload",
            files={"file": ("sample.ab1", b"binary-data", "application/octet-stream")},
        )
        self.assertEqual(resp.status_code, 200)
        path = resp.json()["path"]
        self.assertEqual(
            os.path.dirname(os.path.realpath(path)),
            os.path.realpath(self.uploads_dir),
        )
        self.assertTrue(os.path.exists(path))
        with open(path, "rb") as f:
            self.assertEqual(f.read(), b"binary-data")

    def test_relative_traversal_filename_cannot_escape_uploads_dir(self):
        outside_marker = os.path.join(self.tmp.name, "escaped.txt")
        resp = self.client.post(
            "/upload",
            files={"file": ("../escaped.txt", b"evil", "application/octet-stream")},
        )
        self.assertFalse(
            os.path.exists(outside_marker),
            "traversal filename must never write outside uploads_dir",
        )
        # Either rejected outright, or safely re-based inside uploads_dir.
        if resp.status_code == 200:
            self.assertTrue(
                os.path.realpath(resp.json()["path"]).startswith(
                    os.path.realpath(self.uploads_dir) + os.sep
                )
            )
        else:
            self.assertEqual(resp.status_code, 400)

    def test_absolute_path_filename_cannot_escape_uploads_dir(self):
        target = os.path.join(self.tmp.name, "absolute_escape.txt")
        resp = self.client.post(
            "/upload",
            files={"file": (target, b"evil", "application/octet-stream")},
        )
        self.assertFalse(
            os.path.exists(target),
            "an absolute-path filename must never be honored as-is",
        )
        if resp.status_code == 200:
            self.assertTrue(
                os.path.realpath(resp.json()["path"]).startswith(
                    os.path.realpath(self.uploads_dir) + os.sep
                )
            )
        else:
            self.assertEqual(resp.status_code, 400)


class TestPreviewReadPathTraversal(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)
        self.tmp = tempfile.TemporaryDirectory()
        self.uploads_dir = os.path.join(self.tmp.name, "uploads")
        self.jobs_dir = os.path.join(self.tmp.name, "jobs")
        self.outside_dir = os.path.join(self.tmp.name, "outside")
        for d in (self.uploads_dir, self.jobs_dir, self.outside_dir):
            os.makedirs(d, exist_ok=True)

        self._patchers = [
            patch.object(settings, "uploads_dir", self.uploads_dir),
            patch.object(settings, "jobs_dir", self.jobs_dir),
        ]
        for p in self._patchers:
            p.start()

    def tearDown(self):
        for p in self._patchers:
            p.stop()
        self.tmp.cleanup()

    def test_path_outside_allowed_dirs_is_rejected(self):
        secret = os.path.join(self.outside_dir, "secret.ab1")
        with open(secret, "wb") as f:
            f.write(b"not-really-ab1")

        resp = self.client.get("/preview-read", params={"path": secret})
        self.assertEqual(resp.status_code, 400)

    def test_traversal_out_of_uploads_dir_is_rejected(self):
        secret = os.path.join(self.outside_dir, "secret.ab1")
        with open(secret, "wb") as f:
            f.write(b"not-really-ab1")

        traversal_path = os.path.join(self.uploads_dir, "..", "outside", "secret.ab1")
        resp = self.client.get("/preview-read", params={"path": traversal_path})
        self.assertEqual(resp.status_code, 400)

    def test_nonexistent_path_outside_allowed_dirs_is_still_rejected(self):
        # Containment must be checked before any filesystem existence check,
        # so this doesn't leak "exists vs doesn't" info for arbitrary paths.
        resp = self.client.get("/preview-read", params={"path": "/etc/shadow"})
        self.assertEqual(resp.status_code, 400)

    @patch("services.aligner.get_read_preview")
    def test_path_inside_uploads_dir_is_allowed_through(self, mock_preview):
        mock_preview.return_value = _FAKE_TRACE
        real_file = os.path.join(self.uploads_dir, "sample.ab1")
        with open(real_file, "wb") as f:
            f.write(b"fake-ab1-bytes")

        resp = self.client.get("/preview-read", params={"path": real_file})
        self.assertEqual(resp.status_code, 200)
        mock_preview.assert_called_once()

    @patch("services.aligner.get_read_preview")
    def test_path_inside_jobs_dir_is_allowed_through(self, mock_preview):
        mock_preview.return_value = _FAKE_TRACE
        real_file = os.path.join(self.jobs_dir, "abc123_data", "reads", "sample.ab1")
        os.makedirs(os.path.dirname(real_file), exist_ok=True)
        with open(real_file, "wb") as f:
            f.write(b"fake-ab1-bytes")

        resp = self.client.get("/preview-read", params={"path": real_file})
        self.assertEqual(resp.status_code, 200)
        mock_preview.assert_called_once()


if __name__ == "__main__":
    unittest.main()
