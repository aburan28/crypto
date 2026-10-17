"""Focused controls for sparse-safe lab-browser source discovery."""

import os
import pathlib
import subprocess
import tempfile
import unittest

import build_lab_browser as browser


class RepositoryInputsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name)
        self.git("init", "-q")
        self.git("config", "user.email", "lab-browser-test@example.invalid")
        self.git("config", "user.name", "Lab browser test")
        (self.root / "tracked").mkdir()
        (self.root / "tracked" / "nested").mkdir()
        (self.root / "tracked" / "off.json").write_text('{"source":"head"}\n', encoding="utf-8")
        (self.root / "tracked" / "live.json").write_text('{"source":"old"}\n', encoding="utf-8")
        (self.root / "tracked" / "nested" / "inside.json").write_text(
            '{"source":"nested"}\n', encoding="utf-8"
        )
        self.git("add", "tracked")
        self.git("commit", "-qm", "fixture")
        self.off_oid = self.git("rev-parse", "HEAD:tracked/off.json", capture=True).strip()
        self.git("update-index", "--skip-worktree", "tracked/off.json")
        os.unlink(self.root / "tracked" / "off.json")
        (self.root / "tracked" / "live.json").write_text('{"source":"worktree"}\n', encoding="utf-8")
        (self.root / "tracked" / "untracked.json").write_text(
            '{"source":"untracked"}\n', encoding="utf-8"
        )
        (self.root / "tracked" / ".hidden.json").write_text(
            '{"source":"untracked-hidden"}\n', encoding="utf-8"
        )

    def tearDown(self):
        self.tmp.cleanup()

    def git(self, *args, capture=False):
        result = subprocess.run(
            ["git", *args],
            cwd=self.root,
            check=True,
            stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
            text=True,
        )
        return result.stdout or ""

    def test_sparse_head_and_materialized_worktree_are_overlaid(self):
        source = browser.RepositoryInputs(str(self.root))
        self.assertEqual(source.git_env["GIT_NO_LAZY_FETCH"], "1")
        self.assertEqual(
            source.glob_files("tracked/*.json"),
            ["tracked/live.json", "tracked/off.json"],
        )
        self.assertEqual(source.read_text("tracked/off.json"), '{"source":"head"}\n')
        self.assertEqual(source.read_text("tracked/live.json"), '{"source":"worktree"}\n')
        self.assertEqual(source.size("tracked/off.json"), 18)
        self.assertTrue(source.is_directory("tracked"))
        self.assertEqual(source.child_directories("tracked"), ["nested"])
        self.assertFalse(source.exists("tracked/untracked.json"))
        self.assertTrue(browser._glob_match("tracked/off.json", "tracked/**/*.json"))
        self.assertTrue(
            browser._glob_match("tracked/nested/off.json", "tracked/**/*.json")
        )
        self.assertFalse(
            browser._glob_match("tracked/nested/off.md", "tracked/**/*.json")
        )
        self.assertFalse(
            browser._glob_match("tracked/.hidden.json", "tracked/**/*.json")
        )
        self.assertFalse(
            browser._glob_match("tracked/.hidden/off.json", "tracked/**/*.json")
        )
        self.assertTrue(
            browser._glob_match("tracked/.hidden.json", "tracked/.hidden.json")
        )

    def test_missing_sparse_blob_fails_closed_with_its_path(self):
        source = browser.RepositoryInputs(str(self.root))
        loose = self.root / ".git" / "objects" / self.off_oid[:2] / self.off_oid[2:]
        loose.unlink()
        with self.assertRaisesRegex(OSError, r"HEAD:tracked/off\.json"):
            source.read_text("tracked/off.json")

    def test_non_git_root_fails_closed(self):
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaisesRegex(RuntimeError, "cannot enumerate committed"):
                browser.RepositoryInputs(root)


if __name__ == "__main__":
    unittest.main()
