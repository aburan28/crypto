#!/usr/bin/env python3
"""Offline tests for the ECC2K-130 status snapshot and history merge."""

from __future__ import annotations

import json
import os
import tempfile
import unittest

from render import merge_history
from snapshot import CLAIM_BOUNDARY, FORBIDDEN_PUBLIC_KEYS, assert_public, campaign_state, normalize


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


class SnapshotTests(unittest.TestCase):
    def test_fixture_normalizes_to_public_document(self):
        path = os.path.join(HERE, "testdata", "status.json")
        with open(path, encoding="utf-8") as fh:
            fixture = json.load(fh)
        snapshot = normalize(
            {
                "campaign_id": fixture["campaign_id"],
                "curve_id": fixture["curve_id"],
                "dp_mask_bits": fixture["dp_mask_bits"],
                "campaign_created_at": fixture["campaign_created_at"],
                "dps": fixture["dps"],
                "workers": fixture["workers"],
                "first_dp_at": fixture["first_dp_at"],
                "last_dp_at": fixture["last_dp_at"],
                "dps_last_hour": fixture["dps_last_hour"],
                "dps_last_day": fixture["dps_last_day"],
                "collisions": fixture["collisions"],
                "latest_collision_at": None,
                "per_worker": fixture["per_worker"],
                "hourly": fixture["hourly"],
            },
            "ecc2k-130",
            "fixture",
        )
        self.assertEqual(snapshot["state"], "COLLECTING")
        self.assertEqual(snapshot["claim_boundary"], CLAIM_BOUNDARY)
        self.assertNotIn("point_key", json.dumps(snapshot))
        assert_public(snapshot)

    def test_collision_state(self):
        self.assertEqual(campaign_state({"collisions": 1, "dps": 10, "dps_last_hour": 0}), "COLLISION_RECORDED")
        self.assertEqual(campaign_state({"collisions": 0, "dps": 10, "dps_last_hour": 0}), "IDLE_OR_STALE")
        self.assertEqual(campaign_state({"collisions": 0, "dps": 0, "dps_last_hour": 0}), "EMPTY")

    def test_refuses_walk_secrets(self):
        with self.assertRaises(SystemExit):
            assert_public({"point_key": "00"})
        for key in FORBIDDEN_PUBLIC_KEYS:
            self.assertTrue(key)


class WorkflowTests(unittest.TestCase):
    def test_publish_job_uses_access_keys_and_deploy_key(self):
        path = os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}", text)
        self.assertIn("aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}", text)
        self.assertIn("secrets.RHO_WALKER_SSH_KEY", text)
        self.assertNotIn("role-to-assume", text)

    def test_fetch_script_punches_runner_ip_not_launch_key(self):
        path = os.path.join(HERE, "fetch_via_walker.sh")
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("authorize-security-group-ingress", text)
        self.assertIn("revoke-security-group-ingress", text)
        self.assertIn("checkip.amazonaws.com", text)
        self.assertIn("gha-ecc2k130-status", text)
        self.assertIn("gha_walker.pub", text)
        self.assertNotIn("meow34", text)
        pub = os.path.join(HERE, "gha_walker.pub")
        with open(pub, encoding="utf-8") as fh:
            line = fh.read().strip()
        self.assertTrue(line.startswith("ssh-ed25519 "))
        self.assertIn("ecc2k130-status-gha", line)


class HistoryTests(unittest.TestCase):
    def test_merges_and_dedups(self):
        previous = {
            "points": [
                {"generated_at": "2026-09-12T09:00:00Z", "dps": 100, "dps_last_hour": 10, "collisions": 0, "workers": 1, "state": "COLLECTING"},
                {"generated_at": "2026-09-12T10:00:00Z", "dps": 200, "dps_last_hour": 20, "collisions": 0, "workers": 1, "state": "COLLECTING"},
            ]
        }
        snapshot = {
            "generated_at": "2026-09-12T10:00:00Z",
            "campaign_id": "ecc2k-130",
            "dps": 210,
            "dps_last_hour": 30,
            "collisions": 0,
            "workers": 2,
            "state": "COLLECTING",
        }
        history = merge_history(previous, snapshot, limit=2)
        self.assertEqual(len(history["points"]), 2)
        self.assertEqual(history["points"][-1]["dps"], 210)
        self.assertEqual(history["campaign_id"], "ecc2k-130")

    def test_render_cli(self):
        import render

        with tempfile.TemporaryDirectory() as tmp:
            status = os.path.join(tmp, "status.json")
            history = os.path.join(tmp, "history.json")
            with open(status, "w", encoding="utf-8") as fh:
                json.dump({
                    "generated_at": "2026-09-12T11:00:00Z",
                    "campaign_id": "ecc2k-130",
                    "dps": 1,
                    "dps_last_hour": 1,
                    "collisions": 0,
                    "workers": 1,
                    "state": "COLLECTING",
                }, fh)
            self.assertEqual(render.main(["--status", status, "--history-out", history]), 0)
            with open(history, encoding="utf-8") as fh:
                blob = json.load(fh)
            self.assertEqual(blob["points"][0]["dps"], 1)


if __name__ == "__main__":
    unittest.main()
