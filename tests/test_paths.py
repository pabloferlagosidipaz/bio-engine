# SPDX-License-Identifier: MIT
"""Unit tests for the core.paths safety helpers (safe_filename_join, resolve_within)."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core.paths import PathTraversalError, resolve_within, safe_filename_join


class TestSafeFilenameJoin(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_plain_filename(self):
        result = safe_filename_join(self.base, "sample.ab1")
        self.assertEqual(result, os.path.realpath(os.path.join(self.base, "sample.ab1")))

    def test_relative_traversal_is_stripped_to_basename(self):
        result = safe_filename_join(self.base, "../../etc/passwd")
        self.assertEqual(result, os.path.realpath(os.path.join(self.base, "passwd")))

    def test_absolute_path_is_stripped_to_basename(self):
        result = safe_filename_join(self.base, "/etc/passwd")
        self.assertEqual(result, os.path.realpath(os.path.join(self.base, "passwd")))

    def test_empty_filename_rejected(self):
        with self.assertRaises(PathTraversalError):
            safe_filename_join(self.base, "")

    def test_dot_dot_filename_rejected(self):
        with self.assertRaises(PathTraversalError):
            safe_filename_join(self.base, "..")

    def test_result_status_code_is_400(self):
        try:
            safe_filename_join(self.base, "")
        except PathTraversalError as e:
            self.assertEqual(e.status_code, 400)
        else:
            self.fail("expected PathTraversalError")


class TestResolveWithin(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.allowed = os.path.join(self.tmp.name, "allowed")
        self.outside = os.path.join(self.tmp.name, "outside")
        os.makedirs(self.allowed, exist_ok=True)
        os.makedirs(self.outside, exist_ok=True)

    def tearDown(self):
        self.tmp.cleanup()

    def test_path_inside_allowed_dir_resolves(self):
        target = os.path.join(self.allowed, "sub", "file.ab1")
        result = resolve_within(target, [self.allowed])
        self.assertEqual(result, os.path.realpath(target))

    def test_path_outside_allowed_dirs_rejected(self):
        target = os.path.join(self.outside, "file.ab1")
        with self.assertRaises(PathTraversalError):
            resolve_within(target, [self.allowed])

    def test_traversal_out_of_allowed_dir_rejected(self):
        target = os.path.join(self.allowed, "..", "outside", "file.ab1")
        with self.assertRaises(PathTraversalError):
            resolve_within(target, [self.allowed])

    def test_multiple_allowed_dirs_second_matches(self):
        second_allowed = os.path.join(self.tmp.name, "jobs")
        os.makedirs(second_allowed, exist_ok=True)
        target = os.path.join(second_allowed, "file.ab1")
        result = resolve_within(target, [self.allowed, second_allowed])
        self.assertEqual(result, os.path.realpath(target))

    def test_empty_path_rejected(self):
        with self.assertRaises(PathTraversalError):
            resolve_within("", [self.allowed])

    def test_sibling_directory_with_shared_prefix_is_not_treated_as_contained(self):
        # e.g. "/data/uploads-evil" must not be considered inside "/data/uploads"
        sibling = self.allowed + "-evil"
        os.makedirs(sibling, exist_ok=True)
        target = os.path.join(sibling, "file.ab1")
        with self.assertRaises(PathTraversalError):
            resolve_within(target, [self.allowed])


if __name__ == "__main__":
    unittest.main()
