# SPDX-License-Identifier: MIT
"""
Regression tests for Phase 5: internal exception text (which can contain
filesystem paths or other implementation detail) must be logged server-side
but never forwarded verbatim into the client-facing JSON error response.

Each test forces the underlying operation to fail with a distinctive,
"sensitive-looking" message and asserts that string never appears in the
HTTP response body, while the response still carries a generic, useful
message.
"""
import os
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from fastapi.testclient import TestClient

from core.config import settings
from main import app

_SECRET_DETAIL = "/home/victim/.ssh/id_rsa leaked-internal-detail-12345"


class TestErrorMessageSanitization(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)

    def _assert_generic_and_no_leak(self, resp, expected_status=500):
        self.assertEqual(resp.status_code, expected_status)
        body = resp.json()
        self.assertNotIn(_SECRET_DETAIL, body.get("message", ""))
        self.assertTrue(body.get("message"))  # still a non-empty, useful message

    def test_upload_failure_does_not_leak_exception_text(self):
        with patch.object(settings, "uploads_dir", "/nonexistent/" + _SECRET_DETAIL):
            resp = self.client.post(
                "/upload",
                files={"file": ("sample.ab1", b"data", "application/octet-stream")},
            )
        self._assert_generic_and_no_leak(resp)
        self.assertEqual(resp.json()["message"], "Failed to upload file")

    @patch("services.aligner.get_read_preview")
    def test_preview_read_failure_does_not_leak_exception_text(self, mock_preview):
        mock_preview.side_effect = RuntimeError(_SECRET_DETAIL)
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            uploads_dir = os.path.join(tmp, "uploads")
            os.makedirs(uploads_dir)
            real_file = os.path.join(uploads_dir, "sample.ab1")
            with open(real_file, "wb") as f:
                f.write(b"data")

            with patch.object(settings, "uploads_dir", uploads_dir), \
                 patch.object(settings, "jobs_dir", os.path.join(tmp, "jobs")):
                resp = self.client.get("/preview-read", params={"path": real_file})

        self._assert_generic_and_no_leak(resp)
        self.assertEqual(resp.json()["message"], "Preview failed")

    @patch("utilities.ensembl_hgvs.EnsemblHGVS")
    def test_hgvs_alternatives_failure_does_not_leak_exception_text(self, mock_ensembl_cls):
        mock_ensembl_cls.side_effect = RuntimeError(_SECRET_DETAIL)
        resp = self.client.post(
            "/tools/hgvs/alternatives",
            json={"transcript": "NM_000001", "assembly": "GRCh38", "pos": 1, "ref": "A", "alt": "T"},
        )
        self._assert_generic_and_no_leak(resp)
        self.assertEqual(resp.json()["message"], "HGVS alternatives fetch failed")

    def test_share_job_failure_does_not_leak_exception_text(self):
        with patch("services.job_manager.JobManager.export_job") as mock_export:
            mock_export.side_effect = RuntimeError(_SECRET_DETAIL)
            resp = self.client.post(
                "/jobs/some-job-id/share",
                json={"level": "results_only", "target_folder": None},
            )
        self._assert_generic_and_no_leak(resp)
        self.assertNotIn(_SECRET_DETAIL, resp.json()["message"])

    def test_import_job_failure_does_not_leak_exception_text(self):
        with patch("services.job_manager.JobManager.import_job") as mock_import:
            mock_import.side_effect = RuntimeError(_SECRET_DETAIL)
            resp = self.client.post(
                "/jobs/import",
                json={"source_folder": "/some/folder"},
            )
        self._assert_generic_and_no_leak(resp)
        self.assertEqual(resp.json()["message"], "Failed to import job")


if __name__ == "__main__":
    unittest.main()
