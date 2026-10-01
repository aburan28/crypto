#!/usr/bin/env python3
"""Unit tests for the IC boundary autolab control plane."""
from __future__ import annotations

import json
import os
import stat
import tempfile
import textwrap
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
        self.assertFalse(report["commands"].get("n37_full"))
        self.assertTrue(report["primary_speedup_contract"])
        for row in report["beats"]:
            self.assertFalse(row["primary_speedup_eligible"])
            self.assertTrue(row["primary_speedup_blocker"])


class MeasurementSchemaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.protocol = lab.load_protocol()
        self.ledger = lab.load_ledger(self.protocol)

    def test_vs_rho_incomplete_fails_closed(self) -> None:
        result = lab.validate_claim({"stage": "vs_rho", "n": 37}, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL")
        self.assertIn("timing_class", result["missing_stage_fields"])
        self.assertTrue(result["missing_global_provenance"])

    def test_vs_rho_single_target_online_claim_passes(self) -> None:
        report = {
            "n": 37,
            "candidate_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abc",
            "candidate_manifest_sha256": "123456789abc0000000000000000000000000000000000000000000000000000",
            "workload_id": "abcdef123456",
            "workload_manifest_sha256": "abcdef1234560000000000000000000000000000000000000000000000000000",
            "run_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abcWabcdef123456R1",
            "target_count": 1,
            "ic_target_hash": "target-abc",
            "rho_target_hash": "target-abc",
            "timing_class": "single_target_online_wall",
            "ic_online_wall_ms": 10.0,
            "rho_online_wall_ms": 20.0,
            "online_speedup": 2.0,
            "ic_online_phase_ms": {
                "T_target_query_ms": 1.0,
                "T_target_PDP_ms": 3.0,
                "T_target_relation_check_ms": 1.0,
                "T_target_descent_ms": 4.0,
                "T_target_recovery_check_ms": 1.0,
            },
            "online_interval": {
                "ic_start_event": "first target-dependent IC query after reusable setup",
                "ic_stop_event": "scalar recovered and independently verified",
                "ic_included_stages": ["target_query", "target_PDP", "target_relation_check", "target_descent", "target_recovery_check"],
                "rho_start_event": "first target-dependent rho walk",
                "rho_stop_event": "scalar recovered and independently verified",
                "rho_included_stages": ["walk", "collision", "recovery_check"],
            },
            "same_resource_envelope": True,
            "independent_validation": True,
            "ic_replay_certificate_sha256": "a" * 64,
            "rho_replay_certificate_sha256": "b" * 64,
            "ic_resource_envelope": {"worker_count": 1, "memory_limit_bytes": 1073741824},
            "rho_resource_envelope": {"worker_count": 1, "memory_limit_bytes": 1073741824},
            "ic_scalar_verified": True,
            "rho_scalar_verified": True,
            "rho_policy": {
                "worker_count": 1,
                "walk_policy": "signed_frobenius",
                "collision_policy": "distinguished_point",
                "distinguished_point_memory_bytes": 0,
            },
            "verdict": "DRAFT",
            "claim_boundary": "single public target only",
            "independent_replay_pointer": "research/example",
            "fixture_hash": "target-abc",
            "executable_or_source_hash": {"direct": "def", "rho": "ghi"},
            "host_id": {"node": "test"},
            "resource_caps": {"common_cap_bytes": 1},
            "seeds": {"direct": 1, "rho": 2},
            "claim_boundary_non_claims": ["not key recovery"],
        }
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "PASS", result)

    def test_vs_rho_batch_mismatch_and_legacy_timing_are_rejected(self) -> None:
        import copy

        base = {
            "candidate_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abc",
            "candidate_manifest_sha256": "123456789abc0000000000000000000000000000000000000000000000000000",
            "workload_id": "abcdef123456",
            "workload_manifest_sha256": "abcdef1234560000000000000000000000000000000000000000000000000000",
            "run_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abcWabcdef123456R1",
            "target_count": 1,
            "ic_target_hash": "target-abc",
            "rho_target_hash": "target-abc",
            "timing_class": "single_target_online_wall",
            "ic_online_wall_ms": 10.0,
            "rho_online_wall_ms": 20.0,
            "online_speedup": 2.0,
            "ic_online_phase_ms": {
                "T_target_query_ms": 1.0,
                "T_target_PDP_ms": 3.0,
                "T_target_relation_check_ms": 1.0,
                "T_target_descent_ms": 4.0,
                "T_target_recovery_check_ms": 1.0,
            },
            "online_interval": {
                "ic_start_event": "target computation after reusable setup",
                "ic_stop_event": "verified scalar recovered",
                "ic_included_stages": ["target_query", "target_PDP", "target_relation_check", "target_descent", "target_recovery_check"],
                "rho_start_event": "target walk begins",
                "rho_stop_event": "verified scalar recovered",
                "rho_included_stages": ["walk", "collision", "recovery_check"],
            },
            "same_resource_envelope": True,
            "independent_validation": True,
            "ic_replay_certificate_sha256": "a" * 64,
            "rho_replay_certificate_sha256": "b" * 64,
            "ic_resource_envelope": {"worker_count": 1, "memory_limit_bytes": 1073741824},
            "rho_resource_envelope": {"worker_count": 1, "memory_limit_bytes": 1073741824},
            "ic_scalar_verified": True,
            "rho_scalar_verified": True,
            "rho_policy": {
                "worker_count": 1,
                "walk_policy": "walk",
                "collision_policy": "collision",
                "distinguished_point_memory_bytes": 0,
            },
            "verdict": "DRAFT",
            "claim_boundary": "single public target only",
            "independent_replay_pointer": "research/example",
            "fixture_hash": "target-abc",
            "executable_or_source_hash": {"direct": "def", "rho": "ghi"},
            "host_id": "test",
            "resource_caps": {},
            "seeds": 1,
            "claim_boundary_non_claims": ["not key recovery"],
        }
        variants = []
        batched = copy.deepcopy(base)
        batched["target_count"] = 2
        variants.append(batched)
        for invalid_count in (True, 1.0):
            malformed = copy.deepcopy(base)
            malformed["target_count"] = invalid_count
            variants.append(malformed)
        mismatched = copy.deepcopy(base)
        mismatched["rho_target_hash"] = "another-target"
        variants.append(mismatched)
        legacy = copy.deepcopy(base)
        legacy["timing_class"] = "whole_process_wall"
        variants.append(legacy)
        wrong_ratio = copy.deepcopy(base)
        wrong_ratio["online_speedup"] = 3.0
        variants.append(wrong_ratio)
        unverified = copy.deepcopy(base)
        unverified["rho_scalar_verified"] = False
        variants.append(unverified)
        wrong_phase_total = copy.deepcopy(base)
        wrong_phase_total["ic_online_phase_ms"]["T_target_PDP_ms"] = 4.0
        variants.append(wrong_phase_total)
        wrong_run_id = copy.deepcopy(base)
        wrong_run_id["run_id"] = "legacy-timestamp-run"
        variants.append(wrong_run_id)
        no_independent_validation = copy.deepcopy(base)
        no_independent_validation["independent_validation"] = False
        variants.append(no_independent_validation)
        bad_certificate_hash = copy.deepcopy(base)
        bad_certificate_hash["ic_replay_certificate_sha256"] = "not-a-digest"
        variants.append(bad_certificate_hash)
        mismatched_resources = copy.deepcopy(base)
        mismatched_resources["rho_resource_envelope"]["worker_count"] = 2
        variants.append(mismatched_resources)
        missing_relation_check = copy.deepcopy(base)
        missing_relation_check["online_interval"]["ic_included_stages"].remove("target_relation_check")
        variants.append(missing_relation_check)
        for report in variants:
            result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
            self.assertEqual(result["status"], "FAIL", result)
            self.assertTrue(result["validation_errors"])

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
    def test_matched_targets_are_required_even_when_both_producers_succeed(self):
        row = {'fixture_index':0,'n':13,'a':0,'generator_point_key':['1','2'],
               'subgroup_order':2003,'field_modulus_low_terms':[0,1,3,4],
               'published_q_point_key':['3','4'],'recovered_fixture_scalar':5,'published_fixture_scalar':5}
        direct = [dict(row,kind='point_defined_factor_base'), dict(row,kind='relation_rank_summary',linear_solution_verified=True)]
        rho = [dict(row,kind='rho_public_fixture',verified=True)]
        matched = lab.comparison_integrity(direct,rho,1)
        self.assertEqual(matched['status'],'MATCHED')
        self.assertEqual(len(matched['fixture_hash']),64)
        rho[0]['published_q_point_key'] = [3,6]
        self.assertEqual(lab.comparison_integrity(direct,rho,1)['status'],'INVALID_COMPARISON')
        self.assertIsNone(lab.comparison_integrity(direct,rho,1)['fixture_hash'])

    def test_partial_or_duplicate_workload_fails_closed(self):
        row = {'fixture_index':0,'n':13,'a':0,'generator_point_key':[1,2],
               'subgroup_order':2003,'field_modulus_low_terms':[0,1,3,4],
               'published_q_point_key':[3,4],'recovered_fixture_scalar':5,'published_fixture_scalar':5}
        direct = [dict(row,kind='point_defined_factor_base'), dict(row,kind='relation_rank_summary',linear_solution_verified=True)]
        rho = [dict(row,kind='rho_public_fixture',verified=True)]
        for a,b,count in [(direct,rho,2),(direct+direct,rho,1),([],rho,1)]:
            self.assertEqual(lab.comparison_integrity(a,b,count)['status'],'INVALID_COMPARISON')

    def test_only_single_target_rho_cost_is_comparable(self):
        self.assertIsNone(lab.extract_ic_cost([{'online_charged_ms':1}], 'algorithmic_charged'))
        self.assertIsNone(lab.extract_ic_cost([{'charged_total_ms':1}], 'algorithmic_charged'))
        self.assertEqual(lab.extract_ic_cost([{'full_algorithm_charged_total_ms':9}], 'algorithmic_charged'),9)
        one = {'kind': 'rho_public_fixture', 'total_ms': 2}
        self.assertEqual(lab.extract_rho_cost([one]), 2)
        self.assertIsNone(lab.extract_rho_cost([one, dict(one, fixture_index=1)]))
        self.assertIsNone(lab.extract_rho_cost([{'kind': 'rho_public_fixture', 'total_ms':float('nan')}]))
        self.assertIsNone(lab.extract_ic_cost([{'full_algorithm_charged_total_ms':-1}], 'algorithmic_charged'))

    def test_launch_rejects_more_than_one_target_before_preflight(self):
        args = lab.parser().parse_args([
            'launch', '--beat', 'koblitz.vs_rho.n37_wall', '--fixtures', '2'
        ])
        with self.assertRaisesRegex(lab.AutolabError, 'exactly one target'):
            lab.launch(args)

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
            self.assertEqual(result["status"], "FAIL")
            self.assertIn("target_count", result["missing_stage_fields"])
            self.assertEqual(json.loads(out.read_text())["status"], "FAIL")


class TimingTests(unittest.TestCase):
    """The wall metric must not charge a producer for its first execution.

    A single cold run adds a roughly constant per-process term to both arms,
    which pulls ratios toward parity and so flatters the slower arm.
    """

    def make_producer(self, directory: str, body: str) -> Path:
        path = Path(directory) / "fake_producer"
        path.write_text("#!/usr/bin/env python3\n" + textwrap.dedent(body))
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def test_timing_free_ignores_durations(self) -> None:
        first = json.dumps({"base_hash": "abc", "rank": 3, "setup_ms": 1.5})
        second = json.dumps({"base_hash": "abc", "rank": 3, "setup_ms": 9.9})
        differing = json.dumps({"base_hash": "zzz", "rank": 3, "setup_ms": 1.5})
        self.assertEqual(lab._timing_free(first), lab._timing_free(second))
        self.assertNotEqual(lab._timing_free(first), lab._timing_free(differing))

    def test_warmup_runs_and_is_not_timed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            counter = Path(tmp) / "calls"
            producer = self.make_producer(
                tmp,
                f"""
                import json, sys
                with open({str(counter)!r}, 'a') as handle:
                    handle.write(' '.join(sys.argv[1:]) + '\\n')
                if len(sys.argv) < 2:
                    sys.exit(2)
                print(json.dumps({{'base_hash': 'abc', 'setup_ms': 1.0}}))
                """,
            )
            observed = lab.run_timed(
                [str(producer), "real-arg"],
                env=dict(os.environ),
                cwd=Path(tmp),
                repeats=3,
            )
            calls = counter.read_text().splitlines()
        # One warmup with no arguments, then exactly the timed repeats.
        self.assertEqual(calls[0], "")
        self.assertEqual(calls[1:], ["real-arg"] * 3)
        self.assertEqual(observed["timed_repeats"], 3)
        self.assertEqual(len(observed["whole_process_wall_samples_ms"]), 3)
        self.assertEqual(observed["exit_code"], 0)
        self.assertTrue(observed["stdout_stable"])
        self.assertEqual(
            observed["whole_process_wall_ms"],
            sorted(observed["whole_process_wall_samples_ms"])[1],
        )
        # The warmup's cost is reported, not folded into the result.
        self.assertNotIn(
            observed["cold_start_wall_ms"], observed["whole_process_wall_samples_ms"]
        )

    def test_repeats_stop_once_a_run_exceeds_the_budget(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            producer = self.make_producer(
                tmp,
                """
                import json, sys, time
                if len(sys.argv) < 2:
                    sys.exit(2)
                time.sleep(0.05)
                print(json.dumps({'base_hash': 'abc'}))
                """,
            )
            original = lab.REPEAT_BUDGET_MS
            lab.REPEAT_BUDGET_MS = 1.0
            try:
                observed = lab.run_timed(
                    [str(producer), "real-arg"],
                    env=dict(os.environ),
                    cwd=Path(tmp),
                    repeats=5,
                )
            finally:
                lab.REPEAT_BUDGET_MS = original

        self.assertEqual(observed["timed_repeats"], 1)

    def test_failing_producer_reports_its_exit_code(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            producer = self.make_producer(tmp, "import sys\nsys.exit(3)\n")
            observed = lab.run_timed(
                [str(producer), "real-arg"],
                env=dict(os.environ),
                cwd=Path(tmp),
                repeats=3,
            )

        self.assertEqual(observed["exit_code"], 3)
        self.assertEqual(observed["timed_repeats"], 1)


if __name__ == "__main__":
    unittest.main()
