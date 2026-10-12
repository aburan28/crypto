#!/usr/bin/env python3
"""Tests for scripts/check_ic_leaderboard_paths.py: the glob forms it models, the
workflow text it reads, and that it refuses the forms it does not model."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_ic_leaderboard_paths as chk  # noqa: E402


def match(pat: str, path: str) -> bool:
    return bool(chk.glob_regex(pat).match(path))


class Glob(unittest.TestCase):
    def test_literal_is_exact(self):
        self.assertTrue(match("docs/curves/registry.json", "docs/curves/registry.json"))
        self.assertFalse(match("docs/curves/registry.json", "docs/curves/registry.json.bak"))
        self.assertFalse(match("docs/curves/registry.json", "x/docs/curves/registry.json"))

    def test_dot_is_literal(self):
        self.assertFalse(match("a.json", "aXjson"))

    def test_star_stops_at_slash(self):
        self.assertTrue(match("runs/*/T01.params.json", "runs/k0n41/T01.params.json"))
        self.assertFalse(match("runs/*/T01.params.json", "runs/k0n41/M1/T01.params.json"))
        self.assertTrue(match("tests/ic_*.rs", "tests/ic_framework.rs"))
        self.assertFalse(match("tests/ic_*.rs", "tests/sub/ic_framework.rs"))

    def test_double_star_crosses_slashes(self):
        self.assertTrue(match("docs/ic/**", "docs/ic/runs/a/b.json"))
        self.assertTrue(match("docs/ic/**", "docs/ic/x"))
        self.assertFalse(match("docs/ic/**", "docs/icx/y"))

    def test_double_star_slash_matches_zero_or_more_directories(self):
        pat = "main/**/r*-candidate.price.json"
        self.assertTrue(match(pat, "main/k0n37/M1/r1-candidate.price.json"))
        self.assertTrue(match(pat, "main/r1-candidate.price.json"))
        self.assertFalse(match(pat, "main/k0n37/M1/r1-other.json"))

    def test_unmodelled_forms_are_refused(self):
        for pat in ("a?.json", "a+.json", "[ab].json", "!docs/**"):
            with self.assertRaises(SystemExit):
                chk.glob_regex(pat)


class Workflow(unittest.TestCase):
    def read(self, text: str, event: str) -> list[str]:
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "w.yml"
            f.write_text(text)
            return chk.event_paths(f, event)

    TEXT = """name: x
on:
  pull_request:
    paths:
      - 'a/**'
      - "b/c.json"   # a comment
      - d/e.txt
  push:
    branches: [main]
    paths:
      - 'z/**'
  workflow_dispatch:
jobs:
  j:
    runs-on: ubuntu-latest
"""

    def test_each_event_has_its_own_list(self):
        self.assertEqual(self.read(self.TEXT, "pull_request"), ["a/**", "b/c.json", "d/e.txt"])
        self.assertEqual(self.read(self.TEXT, "push"), ["z/**"])

    def test_absent_event_is_empty(self):
        self.assertEqual(self.read(self.TEXT, "schedule"), [])

    def test_paths_ignore_is_refused(self):
        with self.assertRaises(SystemExit):
            self.read("on:\n  push:\n    paths-ignore:\n      - 'a/**'\n", "push")


if __name__ == "__main__":
    unittest.main()
