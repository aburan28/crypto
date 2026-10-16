#!/usr/bin/env python3
"""Source-only guards for the post-K600 known-answer control."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import reuse_control_supervisor as gate


class ReuseControlSupervisorTests(unittest.TestCase):
    def authorization(self, output: Path) -> dict:
        return {
            "schema": "n83.additional-local-budget/v5",
            "study": gate.STUDY,
            "purpose": "one_known_witness_control_reuse_k600_worker",
            "prior_total_charged_seconds": "17400",
            "prior_k600_audit_sha256": gate.AUDIT_SHA256,
            "reuse_binary_sha256": gate.BINARY_SHA256,
            "additional_authorized_seconds": "1200",
            "output_dir": str(output.resolve()),
            "user_approval_reference": "explicit future approval fixture",
        }

    def test_exact_fresh_grant_is_required(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            good = self.authorization(output)
            self.assertIsNone(gate.authorization_error(good, output))
            for key, value in (
                ("schema", "n83.additional-local-budget/v4"),
                ("prior_total_charged_seconds", "15000"),
                ("prior_k600_audit_sha256", "0" * 64),
                ("reuse_binary_sha256", "0" * 64),
                ("additional_authorized_seconds", "1199"),
                ("output_dir", str(Path(folder) / "different")),
                ("user_approval_reference", ""),
            ):
                self.assertIsNotNone(gate.authorization_error({**good, key: value}, output))

    def test_stale_grant_cannot_launch_or_create_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output, small, envelope = (root / name for name in ("run", "small", "envelope"))
            small.write_bytes(b"s")
            envelope.write_bytes(b"e")
            stale = self.authorization(output)
            stale["prior_total_charged_seconds"] = "15000"
            auth = root / "authorization.json"
            auth.write_text(json.dumps(stale))
            with mock.patch.object(gate.prior, "stage", side_effect=AssertionError("launched")):
                with self.assertRaisesRegex(ValueError, "authorization"):
                    gate.run(auth, output, small, envelope)
            self.assertFalse(output.exists())

    def test_reused_build_and_new_ledger_are_pinned(self):
        self.assertIsNone(gate.prior_evidence_error(gate.ROOT))
        original = gate.prior.sha256_file

        def changed_audit(path: Path) -> str:
            return "0" * 64 if path == gate.ROOT / gate.AUDIT else original(path)

        with mock.patch.object(gate.prior, "sha256_file", side_effect=changed_audit):
            self.assertEqual(gate.prior_evidence_error(gate.ROOT),
                             "pinned K600 evidence bytes differ")
        self.assertEqual(gate.TOTAL_CHARGE, 1200)
        self.assertEqual(gate.PRIOR_TOTAL + gate.TOTAL_CHARGE, 18600)


if __name__ == "__main__":
    unittest.main()
