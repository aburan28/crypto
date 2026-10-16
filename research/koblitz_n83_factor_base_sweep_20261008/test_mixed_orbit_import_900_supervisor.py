"""Fail-closed source controls for the 900-second retained-N83 import gate."""

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import mixed_orbit_import_900_supervisor as gate
import mixed_orbit_supervisor as reference


class ImportGateTests(unittest.TestCase):
    def authorization(self, output: Path) -> dict:
        return {
            "schema": "n83.additional-local-budget/v1", "study": gate.STUDY,
            "purpose": "mixed_orbit_build_import_900_gate",
            "prior_pilot_charged_seconds": "11700",
            "additional_authorized_seconds": "900",
            "output_dir": str(output),
            "user_approval_reference": "synthetic-test-only",
        }

    def test_exact_budget_and_output_binding_precede_docker(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output = root / "out"
            auth = root / "authorization.json"
            record = self.authorization(output)
            record["additional_authorized_seconds"] = "899"
            auth.write_text(json.dumps(record))
            small, envelope = root / "small", root / "envelope"
            small.write_bytes(b"small")
            envelope.write_bytes(b"envelope")
            with patch.object(gate, "image_identity",
                              side_effect=AssertionError("Docker inspected")):
                with self.assertRaisesRegex(ValueError, "additional_authorized_seconds"):
                    gate.gate(auth, output, small, envelope)
            self.assertFalse(output.exists())
            record["additional_authorized_seconds"] = "900"
            self.assertIsNone(gate.authorization_error(record, output))
            self.assertIn("output_dir", gate.authorization_error(record, root / "other"))

    def test_build_and_import_receipts_are_source_bound(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output = root / "out"
            auth = root / "authorization.json"
            auth.write_text(json.dumps(self.authorization(output)))
            small, envelope = root / "small", root / "envelope"
            small.write_bytes(b"small")
            envelope.write_bytes(b"envelope")
            registry = root / "cargo/registry"
            registry.mkdir(parents=True)
            real_hash = gate.sha256_file
            launched = []

            def fake_hash(path):
                if Path(path) == small:
                    return gate.SMALL_SHA256
                if Path(path) == envelope:
                    return gate.ENVELOPE_SHA256
                return real_hash(path)

            def fake_stage(*, name, command, wall, output, container_name):
                launched.append(name)
                self.assertEqual(command[command.index("--network") + 1], "none")
                self.assertEqual(command[command.index("--pull") + 1], "never")
                if name == "build":
                    self.assertEqual(wall, gate.BUILD_WALL)
                    artifact = output / "target/release/examples/koblitz_n83_mixed_orbit_worker"
                    artifact.parent.mkdir(parents=True)
                    artifact.write_bytes(b"\x7fELFsynthetic")
                else:
                    self.assertEqual(wall, gate.VERIFY_WALL)
                    self.assertEqual(command[command.index("--memory") + 1], "4096m")
                    (output / "verify.stdout").write_text(json.dumps({
                        "schema": "n83.mixed-orbit-import/v1", "status": "PASS",
                        "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
                        "small_object_blake3": reference.SMALL_BLAKE3,
                        "envelope_object_blake3": reference.ENVELOPE_BLAKE3,
                        "small_points": 196212, "small_signed_orbits": 1182,
                        "residual_points": 143756, "residual_signed_orbits": 866,
                        "solver_search_executed": False,
                    }))
                return {"name": name, "exit_code": 0, "timed_out": False,
                        "cleanup_verified": True, "launch_error": None,
                        "container_state": {"ExitCode": 0, "OOMKilled": False}}

            with patch.dict(os.environ, {"CARGO_HOME": str(root / "cargo")}), \
                 patch.object(gate, "clean_commit", return_value="a" * 40), \
                 patch.object(gate, "sha256_file", side_effect=fake_hash), \
                 patch.object(gate, "image_identity", return_value=gate.IMAGE), \
                 patch.object(gate, "stage", side_effect=fake_stage):
                result = gate.gate(auth, output, small, envelope)
            self.assertEqual(launched, ["build", "verify"])
            self.assertEqual(result["status"], "PASS_import_only")
            self.assertEqual(result["conservative_charge_seconds"], 900)
            self.assertIsNone(result["total_index_calculus_runtime_ms"])
            self.assertEqual(json.loads((output / "budget.json").read_text())
                             ["remaining_additional_seconds"], "0")


if __name__ == "__main__":
    unittest.main()
