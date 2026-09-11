#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import koblitz_stage26_affinity_inputs as inputs
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_stage26_affinity_cell as cell


class Stage26Tests(unittest.TestCase):
    def test_truth_keys_are_rejected_recursively(self) -> None:
        inputs.validate_no_truth({"factor_base": {"algebraic": True}})
        with self.assertRaises(inputs.Stage26InputError):
            inputs.validate_no_truth({"nested": [{"known_witness": [1, 2, 3]}]})

    def test_cell_seal_verifies_and_tamper_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = {
                "schema": cell.SCHEMA,
                "status": "complete_truth_free_single_cpu_cell",
                "cell_id": "n31-l5-m3-standard-a1-f0",
                "instances": 40,
                "backend_outcomes": 120,
                "truth_labels_present": False,
                "known_witnesses_present": False,
                "koblitz_index_calculus_sota": False,
            }
            phase_b.write_json_new(root / "result.json", result)
            inventory = phase_b.all_regular_inventory(root, {"result-seal.json"})
            payload = {
                "schema": cell.SEAL_SCHEMA,
                "status": "cell_frozen",
                "cell_id": result["cell_id"],
                "packet_inventory_sha256": "c" * 64,
                "inventory": inventory,
                "inventory_sha256": phase_b.canonical_sha256(inventory),
            }
            phase_b.write_json_new(
                root / "result-seal.json",
                {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)},
            )
            self.assertEqual(cell.verify(root)["status"], "verified")
            result["instances"] = 39
            (root / "result.json").write_text(json.dumps(result))
            with self.assertRaises(cell.Stage26Error):
                cell.verify(root)

    @unittest.skipUnless(sys.platform.startswith("linux"), "Linux /proc control")
    def test_process_tree_sampler_sees_current_process(self) -> None:
        self.assertGreater(cell.proc_tree_rss(cell.os.getpid()), 0)


if __name__ == "__main__":
    unittest.main()
