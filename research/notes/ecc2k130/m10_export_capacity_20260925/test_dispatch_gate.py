"""Fail-closed controls for the hosted one-shot capacity admission."""
from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import dispatch_gate


class DispatchGateTests(unittest.TestCase):
    def test_latest_review_dismissal_supersedes_approval(self):
        rows = [
            {"user": {"login": "reviewer"}, "state": "APPROVED",
             "submitted_at": "2026-09-29T01:00:00Z"},
            {"user": {"login": "reviewer"}, "state": "DISMISSED",
             "submitted_at": "2026-09-29T02:00:00Z"},
        ]
        self.assertEqual(dispatch_gate.latest_reviews(rows)["reviewer"]["state"],
                         "DISMISSED")

    def test_prior_capacity_job_blocks_second_attempt(self):
        def api(path, _token):
            if path.endswith("/jobs?per_page=100"):
                return {"total_count": 1, "jobs": [
                    {"id": 99, "name": "capacity", "conclusion": "failure"}]}
            return {"total_count": 1, "workflow_runs": [
                {"id": 12, "head_branch": dispatch_gate.BRANCH}]}
        with patch.object(dispatch_gate, "api", side_effect=api):
            with self.assertRaisesRegex(RuntimeError, "prior capacity job exists"):
                dispatch_gate.check_prior_capacity("owner/repo", dispatch_gate.BRANCH, 13, "token")

    def test_missing_auth_archives_refusal_without_child(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "gate"
            with patch.dict(os.environ, {}, clear=True):
                with self.assertRaises(KeyError):
                    dispatch_gate.admit(out)
            receipt = json.loads((out / "predispatch.json").read_text())
            self.assertEqual(receipt["decision"], "REFUSED")
            self.assertEqual(receipt["measurement_children"], 0)
            self.assertEqual(receipt["error_type"], "KeyError")


    def test_follow_on_pr_requires_independent_exact_head_review(self):
        head = "a" * 40
        pr_number = 936
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            event_path = folder / "event.json"
            event_path.write_text(json.dumps({
                "action": "labeled",
                "label": {"name": dispatch_gate.LABEL},
                "number": pr_number,
                "pull_request": {
                    "number": pr_number,
                    "head": {"ref": dispatch_gate.BRANCH,
                             "repo": {"full_name": "owner/repo"}, "sha": head},
                },
            }))
            environment = {
                "GITHUB_TOKEN": "test-token", "GITHUB_REPOSITORY": "owner/repo",
                "GITHUB_RUN_ID": "100", "GITHUB_RUN_ATTEMPT": "1",
                "GITHUB_EVENT_PATH": str(event_path),
            }
            live = {"number": pr_number, "state": "open", "draft": False,
                    "head": {"sha": head, "ref": dispatch_gate.BRANCH},
                    "user": {"id": 1}}
            def api(path, _token):
                if path.endswith("/reviews?per_page=100"):
                    return []
                return live
            out = folder / "gate"
            with patch.dict(os.environ, environment, clear=True), \
                 patch.object(dispatch_gate, "api", side_effect=api), \
                 patch.object(dispatch_gate.subprocess, "check_output", return_value=head + "\n"):
                with self.assertRaisesRegex(AssertionError, "independent exact-head approval is missing"):
                    dispatch_gate.admit(out)
            receipt = json.loads((out / "predispatch.json").read_text())
            self.assertEqual(receipt["decision"], "REFUSED")
            self.assertEqual(receipt["measurement_children"], 0)


if __name__ == "__main__":
    unittest.main()
