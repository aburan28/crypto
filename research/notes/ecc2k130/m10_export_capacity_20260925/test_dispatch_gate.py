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


if __name__ == "__main__":
    unittest.main()
