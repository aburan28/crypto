#!/usr/bin/env python3
"""Test fail-closed interpretation of partial stage receipts."""

import json
import unittest

from check_trace import MARKER_STAGES, classify, selftest_receipt


def marker(index):
    return {
        "schema": "n83.fixed-witness-stage-marker/v1",
        "sequence": index,
        "stage": MARKER_STAGES[index],
        "elapsed_process_wall_ms_L0": float(index),
        "input_blake3": None if index == 0 else "a" * 64,
        "solver_size": None if index < 2 else {
            "variables": 100 + index, "clauses": 200 + index,
            "xors": 10 + index,
        },
    }


def trace(rows):
    return ("\n".join(json.dumps(row) for row in rows) + "\n").encode()


class StageTraceTests(unittest.TestCase):
    def test_selftest_transport_receipt(self):
        row = selftest_receipt()
        self.assertEqual(row["status"], "PASS_marker_transport_only")
        self.assertFalse(row["solver_search_executed"])

    def test_timeout_before_solver_marker_is_censored(self):
        row = classify(b"", trace([marker(i) for i in range(6)]), 124)
        self.assertEqual(row["status"], "CENSORED_before_solver_entry_marker")
        self.assertEqual(row["last_completed_stage"], "units_pinned")
        self.assertFalse(row["solver_entry_marker_seen"])
        self.assertIsNone(row["natural_relations"])

    def test_timeout_after_solver_marker_is_still_censored(self):
        row = classify(b"", trace([marker(i) for i in range(7)]), 124)
        self.assertEqual(row["status"], "CENSORED_after_solver_entry_marker")
        self.assertTrue(row["solver_entry_marker_seen"])
        self.assertIsNone(row["selected_best_total_runtime"])

    def test_broken_sequence_and_truncated_json_fail_closed(self):
        rows = [marker(i) for i in range(3)]
        rows[2]["sequence"] = 3
        with self.assertRaisesRegex(ValueError, "stage order"):
            classify(b"", trace(rows), 124)
        with self.assertRaisesRegex(ValueError, "malformed stage JSON"):
            classify(b"", trace([marker(0)]) + b'{"schema":', 124)

    def test_input_or_clock_change_fails_closed(self):
        rows = [marker(i) for i in range(3)]
        rows[2]["input_blake3"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "input digest changed"):
            classify(b"", trace(rows), 124)
        rows[2]["input_blake3"] = "a" * 64
        rows[2]["elapsed_process_wall_ms_L0"] = 0.5
        with self.assertRaisesRegex(ValueError, "stage wall time"):
            classify(b"", trace(rows), 124)

    def test_unsat_receipt_is_encoding_diagnostic(self):
        result = {"schema": "n83.fixed-witness-sat-control/v2",
                  "stage_markers_emitted": 9, "input_blake3": "a" * 64,
                  "status": "UNSAT_fixed_witness_encoding"}
        row = classify(json.dumps(result).encode(),
                       trace([marker(i) for i in range(9)]), 0)
        self.assertEqual(row["status"], "ENCODING_CONTRADICTION_requires_diagnosis")
        self.assertIsNone(row["natural_relations"])


if __name__ == "__main__":
    unittest.main()
