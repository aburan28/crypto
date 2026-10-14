#!/usr/bin/env python3
"""Fail-closed metadata checks for frozen labels and active added text."""
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import check_curve_names as names


REPORT = "research/f4_f5_history_20261012/reports/run-20261012T0009Z/REPORT.md"


class ArchivedLabelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = (names.REPO / REPORT).read_text()
        cls.resolver = names.Resolver()

    def test_exact_frozen_labels_resolve_to_the_reviewed_models(self):
        hits = names.find_retired(self.text, REPORT, self.resolver)
        self.assertEqual(len(hits), 2)
        for _, _, found, replacement in hits:
            self.assertTrue(names.archived_label(REPORT, self.text, found, replacement))
        with patch.object(names, "tracked", return_value=[REPORT]):
            self.assertEqual(names.lint(self.resolver, {}, False), [])

    def test_changed_hash_unknown_label_wrong_model_and_different_path_are_rejected(self):
        label = "K_1 / 2^23"
        model = "icv1-f2m23-tm5197-1f85e9e1"
        self.assertFalse(names.archived_label(REPORT, self.text + "\n", label, model))
        self.assertFalse(names.archived_label(REPORT, self.text, "K_1 / 2^19", model))
        self.assertFalse(names.archived_label(REPORT, self.text, label, "icv1-f2m17-tm101-00378d4e"))
        self.assertFalse(names.archived_label("ordinary/README.md", self.text, label, model))
        self.assertFalse(names.archived_label(REPORT, self.text, label, None))

    def test_diff_keeps_only_the_exact_archived_labels(self):
        diff = (f"+++ b/{REPORT}\n@@ -0,0 +1,3 @@\n"
                "+K_1 / 2^23\n+K_1 / 2^17\n+K_1 / 2^19\n")
        with patch.object(names, "git", side_effect=["", diff]):
            problems = names.lint_diff("base", self.resolver)
        self.assertEqual(len(problems), 1)
        self.assertIn("K_1 / 2^19", problems[0])

    def diff_problems(self, rel, text, pins=None, new_file=""):
        lines = "".join("+" + line + "\n" for line in text.splitlines())
        diff = f"+++ b/{rel}\n@@ -0,0 +1,{len(text.splitlines())} @@\n{lines}"
        with patch.object(names, "git", side_effect=[new_file, diff]):
            return names.lint_diff("base", self.resolver, pins)

    def test_diff_preserves_each_existing_frozen_directory(self):
        for directory in names.FROZEN_DIRS:
            with self.subTest(directory=directory):
                rel = f"research/study/{directory}/snapshot.rs"
                self.assertEqual(self.diff_problems(rel, "K_1 / 2^23"), [])

    def test_diff_still_rejects_active_text_and_nearby_directory_names(self):
        for rel in ("research/study/live.rs", "research/study/archives_notes/live.rs",
                    "research/study/checks/live.rs", "research/study/archives.rs",
                    "docs/current.md"):
            with self.subTest(rel=rel):
                problems = self.diff_problems(rel, "K_1 / 2^23")
                self.assertEqual(len(problems), 1)
                self.assertIn("added line uses retired curve name", problems[0])

    def test_diff_hash_exemption_requires_exact_contents_and_local_manifest(self):
        rel = "research/study/snapshot.rs"
        text = "K_1 / 2^23\n"
        digest = hashlib.sha256(text.encode()).hexdigest()
        pins = {digest: {"research/study"}}
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            path = repo / rel
            path.parent.mkdir(parents=True)
            path.write_text(text)
            with patch.object(names, "REPO", repo):
                self.assertEqual(self.diff_problems(rel, text, pins), [])
                unrelated_pins = {digest: {"research/other-study"}}
                self.assertEqual(len(self.diff_problems(rel, text, unrelated_pins)), 1)
                path.write_text(text + "K_1 / 2^19\n")
                self.assertEqual(len(self.diff_problems(rel, path.read_text(), pins)), 2)

    def test_diff_keeps_retired_filenames_rejected_outside_frozen_evidence(self):
        archive = "research/study/archives/k1n23.rs"
        active = "research/study/k1n23.rs"
        self.assertEqual(self.diff_problems(archive, "", new_file=archive), [])
        problems = self.diff_problems(active, "", new_file=active)
        self.assertEqual(len(problems), 1)
        self.assertIn("new file named by a retired curve stem", problems[0])


if __name__ == "__main__":
    unittest.main()
