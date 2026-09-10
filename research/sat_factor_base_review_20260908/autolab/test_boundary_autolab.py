#!/usr/bin/env python3
"""Unit tests for the IC boundary autolab control plane."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import boundary_autolab as lab


class LedgerContractTests(unittest.TestCase):
    def test_protocol_and_ledger_load(self) -> None:
        protocol = lab.load_protocol()
        ledger = lab.load_ledger(protocol)
        self.assertEqual(protocol["task_id"], lab.TASK_ID)
        self.assertEqual(ledger["schema_version"], 2)
        self.assertTrue(ledger["measurement_schema"]["fail_closed"])
        self.assertIn("koblitz.vs_rho.n37_wall", protocol["beats"])
        self.assertIn("koblitz.vs_rho.n41_charged", protocol["beats"])

    def test_plan_lists_priority_beats(self) -> None:
        protocol = lab.load_protocol()
        ledger = lab.load_ledger(protocol)
        report = lab.plan(protocol, ledger)
        beat_ids = [row["beat_id"] for row in report["beats"]]
        self.assertIn("smoke.koblitz.vs_rho.n13", beat_ids)
        self.assertIn("koblitz.vs_rho.n37_wall", beat_ids)
        self.assertTrue(report["agent_priorities"])


class MeasurementSchemaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.protocol = lab.load_protocol()
        self.ledger = lab.load_ledger(self.protocol)

    def test_vs_rho_incomplete_fails_closed(self) -> None:
        result = lab.validate_claim({"stage": "vs_rho", "n": 37}, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("timing_class", result["missing_stage_fields"])
        self.assertTrue(result["missing_global_provenance"])

    def test_vs_rho_complete_passes(self) -> None:
        report = {
            "n": 37,
            "timing_class": "whole_process_wall",
            "ic_cost": 1.0,
            "rho_cost": 2.0,
            "automorphism_discount": {"A": 74, "formula": "sqrt(2*n)"},
            "all_stages_charged_same_series": True,
            "verdict": "DRAFT",
            "claim_boundary": "synthetic only",
            "independent_replay_pointer": "research/example",
            "fixture_hash": "abc",
            "executable_or_source_hash": "def",
            "host_id": {"node": "test"},
            "resource_caps": {"common_cap_bytes": 1},
            "seeds": {"direct": 1},
            "claim_boundary_non_claims": ["not key recovery"],
        }
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "PASS", result)

    def test_decomposition_requires_ffd(self) -> None:
        report = {
            "n": 31,
            "m": 2,
            "dim": 16,
            "unknowns": 32,
            "system_degree": "quadratic",
            "eq_var_ratio": 1.5,
            "oracle_class": "sat",
            "median_ms_per_target": 10.0,
            "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "budget": "1h"},
            "fixture_hash": "a",
            "executable_or_source_hash": "b",
            "host_id": "h",
            "resource_caps": {},
            "seeds": 1,
            "claim_boundary_non_claims": ["synthetic"],
        }
        result = lab.validate_claim(report, stage="decomposition", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("ffd_or_degree_of_regularity", result["missing_stage_fields"])


class HelperTests(unittest.TestCase):
    def test_seed_is_deterministic(self) -> None:
        self.assertEqual(
            lab.seed_for("koblitz.vs_rho.n37_wall", "direct", 0),
            lab.seed_for("koblitz.vs_rho.n37_wall", "direct", 0),
        )
        self.assertNotEqual(
            lab.seed_for("koblitz.vs_rho.n37_wall", "direct", 0),
            lab.seed_for("koblitz.vs_rho.n37_wall", "rho", 0),
        )

    def test_claim_check_cli_roundtrip(self) -> None:
        report = {
            "stage": "vs_rho",
            "n": 13,
            "timing_class": "algorithmic_charged",
            "ic_cost": 10,
            "rho_cost": 20,
            "automorphism_discount": "sqrt(2n)",
            "all_stages_charged_same_series": True,
            "verdict": "TEST",
            "claim_boundary": "synthetic",
            "independent_replay_pointer": "x",
            "fixture_hash": "a",
            "executable_hash": "b",
            "host_id": "h",
            "resource_caps": {},
            "seed": 1,
            "non_claims": ["not key recovery"],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.json"
            out = Path(tmp) / "check.json"
            path.write_text(json.dumps(report))
            args = lab.parser().parse_args(
                ["claim-check", "--report", str(path), "--stage", "vs_rho", "--out", str(out)]
            )
            result = lab.claim_check(args)
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(json.loads(out.read_text())["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
