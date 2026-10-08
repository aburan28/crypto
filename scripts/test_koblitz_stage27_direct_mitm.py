#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_stage27_direct_mitm as stage27


class Stage27Tests(unittest.TestCase):
    def process(self) -> dict:
        return {
            "returncode": 0,
            "timed_out": False,
            "metrics": {
                "wall_seconds": 1.0,
                "user_seconds": 0.8,
                "system_seconds": 0.1,
                "total_core_seconds": 0.9,
                "single_core_seconds": 0.9,
                "peak_rss_bytes": 1024,
                "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
            },
        }

    def report(self, status: str) -> dict:
        return {
            "schema": "koblitz_pdp_isolated_backend.v1",
            "backend": "direct-mitm",
            "status": status,
            "source_instance_verified": True,
            "regenerated_source_exact": True,
            "exhaustive": True,
            "factor_points": 64,
            "pair_entries": 2000,
            "group_additions": 2080,
            "witness_indices": [1, 2, 3] if status == "sat" else None,
        }

    def test_sat_and_unsat_terminals(self) -> None:
        for status in ("sat", "unsat"):
            stage27.validate_terminal(
                {
                    "solver": "direct-mitm",
                    "status": status,
                    "source_witness_valid": True if status == "sat" else None,
                    "process": self.process(),
                    "backend_report": self.report(status),
                },
                "blind-id",
            )

    def test_invalid_sat_witness_fails(self) -> None:
        row = {
            "solver": "direct-mitm",
            "status": "sat",
            "source_witness_valid": False,
            "process": self.process(),
            "backend_report": self.report("sat"),
        }
        with self.assertRaises(stage27.Stage27Error):
            stage27.validate_terminal(row, "blind-id")

    def test_seal_verifier_detects_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = {
                "schema": stage27.SCHEMA,
                "status": "complete_truth_free_single_cpu_direct_mitm_cell",
                "cell_id": "n31-l5-m3-standard-a1-f0",
                "instances": 40,
                "outcomes": 40,
                "status_counts": {"sat": 20, "unsat": 20},
                "truth_labels_present": False,
                "known_witnesses_present": False,
                "koblitz_index_calculus_sota": False,
            }
            phase_b.write_json_new(root / "result.json", result)
            inventory = phase_b.all_regular_inventory(root, {"result-seal.json"})
            payload = {
                "schema": stage27.SEAL_SCHEMA,
                "status": "direct_mitm_cell_frozen",
                "cell_id": result["cell_id"],
                "packet_inventory_sha256": stage27.PACKET_SHA256,
                "inventory": inventory,
                "inventory_sha256": phase_b.canonical_sha256(inventory),
            }
            phase_b.write_json_new(root / "result-seal.json", {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)})
            self.assertEqual(stage27.verify(root)["status"], "verified")
            (root / "result.json").write_text("{}\n")
            with self.assertRaises(stage27.Stage27Error):
                stage27.verify(root)


if __name__ == "__main__":
    unittest.main()
