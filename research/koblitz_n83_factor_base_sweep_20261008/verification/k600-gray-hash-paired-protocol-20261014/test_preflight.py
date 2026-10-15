"""Controls for the source-only K600 paired gate protocol."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import preflight


class PairedProtocolTest(unittest.TestCase):
    def test_changed_parent_audit_is_rejected_before_planning(self):
        source, expected = preflight.PINNED_SHA256["live_audit"]
        with TemporaryDirectory() as directory:
            changed = Path(directory) / "audit.json"
            value = json.loads(source.read_text())
            value["pairs"]["gray_pool"]["residual_signed_orbits"] += 1
            changed.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError, "pinned bytes changed"):
                preflight.checked_json(changed, expected)

    def test_saved_arms_share_target_base_and_resource_contract(self):
        value = json.loads(preflight.PLAN.read_text())
        self.assertEqual(value["status"], "SOURCE_READY_AUTHORIZATION_PENDING")
        self.assertFalse(value["new_n83_solver_search_executed"])
        self.assertIsNone(value["selected_best_total_runtime"])
        self.assertEqual(value["additional_seconds_requested_not_authorized"], 4800)
        gray, hashed = value["arms"]
        self.assertEqual((gray["residual_signed_orbits"], hashed["residual_signed_orbits"]),
                         (97, 600))
        self.assertEqual(gray["base_object_blake3"], hashed["base_object_blake3"])
        for arm in (gray, hashed):
            self.assertEqual(arm["worker_search_argv"][1], "run-cross-pilot")
            self.assertEqual(arm["worker_search_argv"][2], preflight.BASE)
            self.assertEqual(arm["worker_search_argv"][-5:-1], ["5000", "32", "3", "2"])


if __name__ == "__main__":
    unittest.main()
