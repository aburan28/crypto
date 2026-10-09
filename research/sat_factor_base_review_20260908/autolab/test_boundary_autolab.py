#!/usr/bin/env python3
"""Unit tests for the IC boundary autolab control plane."""
from __future__ import annotations

import json
import io
import os
import stat
import tempfile
import textwrap
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

import boundary_autolab as lab
import op_accounting


def valid_operation_accounting(ic_ops: float = 100.0, rho_ops: float = 400.0) -> dict:
    return {
        "unit_assumption": "IC probes and rho steps are distinct native counters",
        "comparison_status": "native_counters_only",
        "operation_units": {"ic": "target probes", "rho": "walk steps"},
        "ic_online_operations": ic_ops,
        "rho_online_operations": rho_ops,
        "rho_per_ic_native_counter": rho_ops / ic_ops,
        "ops_speedup_online": None,
    }


class LedgerContractTests(unittest.TestCase):
    def test_protocol_and_ledger_load(self) -> None:
        protocol = lab.load_protocol()
        ledger = lab.load_ledger(protocol)
        self.assertEqual(protocol["task_id"], lab.TASK_ID)
        self.assertEqual(ledger["schema_version"], 2)
        self.assertTrue(ledger["measurement_schema"]["fail_closed"])
        self.assertIn("koblitz.vs_rho.n37_wall", protocol["beats"])
        self.assertIn("koblitz.vs_rho.n41_charged", protocol["beats"])
        self.assertIn("koblitz.factor_base.n53", protocol["beats"])
        self.assertNotIn("rho", protocol["beats"]["koblitz.factor_base.n53"])
        self.assertEqual(
            protocol["beats"]["koblitz.factor_base.n53"]["stage"], "factor_base"
        )

    def test_plan_lists_priority_beats(self) -> None:
        protocol = lab.load_protocol()
        ledger = lab.load_ledger(protocol)
        report = lab.plan(protocol, ledger)
        beat_ids = [row["beat_id"] for row in report["beats"]]
        self.assertIn("smoke.koblitz.vs_rho.n13", beat_ids)
        self.assertIn("koblitz.vs_rho.n37_wall", beat_ids)
        self.assertTrue(report["agent_priorities"])
        self.assertNotIn("n37_full", report["commands"])
        self.assertIn("n37_single_target", report["commands"])

    def test_status_overlays_rejected_paired_target_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            (run / "artifacts").mkdir()
            original_state = {
                "status": "PENDING_INDEPENDENT_VALIDATION",
                "claim_check": "PASS",
            }
            (run / "state.json").write_text(json.dumps(original_state))
            (run / "artifacts/paired_target_audit.json").write_text(json.dumps({
                "status": "PAIRING_REJECTED",
                "speedup_claim_valid": False,
                "reason": "IC and rho public points differ",
            }))
            output = io.StringIO()
            with mock.patch.object(lab, "resolve_run", return_value=run):
                with redirect_stdout(output):
                    lab.status(mock.Mock(run_id="fixture"))

            effective = json.loads(output.getvalue())
            self.assertEqual(effective["status"], "PAIRING_FAILURE")
            self.assertEqual(effective["claim_check"], "FAIL")
            self.assertEqual(effective["recorded_status"], "PENDING_INDEPENDENT_VALIDATION")
            self.assertEqual(effective["recorded_claim_check"], "PASS")
            self.assertEqual(effective["pairing_audit_status"], "PAIRING_REJECTED")
            self.assertEqual(
                json.loads((run / "state.json").read_text()), original_state
            )


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
            "n_or_bits": 37,
            "timing_class": "single_target_online",
            "ic_cost": 10.0,
            "rho_cost": 20.0,
            "automorphism_discount": {"A": 74, "formula": "sqrt(2*n)"},
            "all_stages_charged_same_series": True,
            "candidate_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abc",
            "candidate_manifest_sha256": "123456789abc0000000000000000000000000000000000000000000000000000",
            "workload_id": "abcdef123456",
            "workload_manifest_sha256": "abcdef1234560000000000000000000000000000000000000000000000000000",
            "run_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abcWabcdef123456R1",
            "target_count": 1,
            "ic_target_hash": "target-abc",
            "rho_target_hash": "target-abc",
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
            "paired_target": {
                "ic_public_q": [17, 23],
                "rho_public_q": [17, 23],
                "same_public_point": True,
            },
            "ic_verified": True,
            "rho_verified": True,
            "ic_online_ms": 10.0,
            "rho_online_ms": 20.0,
            "rho_online_phase_ms": {"target_walk": 20.0},
            "ic_online_interval": "target work only",
            "rho_online_interval": "target walk through replay",
            "fixture_hash": "target-abc",
            "executable_or_source_hash": {"direct": "def", "rho": "ghi"},
            "host_id": {"node": "test"},
            "resource_caps": {"common_cap_bytes": 1},
            "seeds": {"direct": 1, "rho": 2},
            "claim_boundary_non_claims": ["not key recovery"],
            "operation_accounting": valid_operation_accounting(),
        }
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "PASS", result)

        del report["operation_accounting"]
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "PASS", result)

        report["operation_accounting"] = valid_operation_accounting()
        report["operation_accounting"]["ops_speedup_online"] = 40.0
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL")
        self.assertIn(
            "operation_accounting.ops_speedup_online requires a calibrated common unit and receipt",
            result["pairing_errors"],
        )

        report["operation_accounting"] = valid_operation_accounting()
        del report["operation_accounting"]["operation_units"]
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL")
        self.assertIn(
            "operation_accounting.operation_units must name the ic and rho units",
            result["pairing_errors"],
        )

        report["operation_accounting"] = valid_operation_accounting()
        report["ic_online_phase_ms"]["T_target_query_ms"] = "invalid"
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL")
        self.assertIn(
            "ic_online_phase_ms must contain finite nonnegative costs",
            result["pairing_errors"],
        )

    def test_operation_accounting_rejects_nonfinite_counts_and_ratios(self) -> None:
        for key, value in (
            ("ic_online_operations", float("nan")),
            ("rho_online_operations", float("inf")),
            ("rho_per_ic_native_counter", float("nan")),
        ):
            with self.subTest(key=key):
                accounting = valid_operation_accounting()
                accounting[key] = value
                self.assertTrue(lab.operation_accounting_errors(accounting))

        accounting = valid_operation_accounting()
        accounting["comparison_status"] = "calibrated_common_unit"
        accounting["operation_units"] = {"ic": "group additions", "rho": "group additions"}
        accounting["calibration_receipt"] = "calibration.json"
        accounting["ops_speedup_online"] = float("inf")
        self.assertTrue(lab.operation_accounting_errors(accounting))

    def test_vs_rho_different_public_points_fail(self) -> None:
        report = {
            "n_or_bits": 41,
            "timing_class": "single_target_online",
            "ic_cost": 10.0,
            "rho_cost": 20.0,
            "automorphism_discount": {"A": 82, "formula": "sqrt(2*n)"},
            "all_stages_charged_same_series": True,
            "verdict": "DRAFT",
            "claim_boundary": "synthetic only",
            "independent_replay_pointer": "research/example",
            "target_count": 1,
            "paired_target": {
                "ic_public_q": [17, 23],
                "rho_public_q": [19, 29],
                "same_public_point": False,
            },
            "ic_verified": True,
            "rho_verified": True,
            "ic_online_ms": 10.0,
            "rho_online_ms": 20.0,
            "online_speedup": 2.0,
            "ic_online_phase_ms": {"target_pdp": 10.0},
            "rho_online_phase_ms": {"target_walk": 20.0},
            "ic_online_interval": "target work only",
            "rho_online_interval": "target walk through replay",
            "fixture_hash": "abc",
            "executable_or_source_hash": "def",
            "host_id": {"node": "test"},
            "resource_caps": {"common_cap_bytes": 1},
            "seeds": {"direct": 1, "rho": 2},
            "claim_boundary_non_claims": ["not key recovery"],
            "operation_accounting": valid_operation_accounting(),
        }
        result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(result["status"], "FAIL", result)
        self.assertIn("IC and rho public targets differ", result["pairing_errors"])
    def test_vs_rho_batch_mismatch_and_legacy_timing_are_rejected(self) -> None:
        import copy

        base = {
            "n_or_bits": 37,
            "ic_cost": 10.0,
            "rho_cost": 20.0,
            "automorphism_discount": {"A": 74, "formula": "sqrt(2*n)"},
            "all_stages_charged_same_series": True,
            "candidate_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abc",
            "candidate_manifest_sha256": "123456789abc0000000000000000000000000000000000000000000000000000",
            "workload_id": "abcdef123456",
            "workload_manifest_sha256": "abcdef1234560000000000000000000000000000000000000000000000000000",
            "run_id": "IC1N37Ckb1fb64PDP5f4RCwalkLAbwTDdirectISO0h123456789abcWabcdef123456R1",
            "target_count": 1,
            "ic_target_hash": "target-abc",
            "rho_target_hash": "target-abc",
            "timing_class": "single_target_online",
            "ic_online_wall_ms": 10.0,
            "rho_online_wall_ms": 20.0,
            "ic_online_ms": 10.0,
            "rho_online_ms": 20.0,
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
            "ic_verified": True,
            "rho_verified": True,
            "paired_target": {
                "ic_public_q": [17, 23],
                "rho_public_q": [17, 23],
                "same_public_point": True,
            },
            "rho_online_phase_ms": {"target_walk": 20.0},
            "ic_online_interval": "first target query through independent replay",
            "rho_online_interval": "first target walk through independent replay",
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
        self.assertEqual(lab.validate_claim(base, stage="vs_rho", ledger=self.ledger)["status"], "PASS")
        for index, report in enumerate(variants):
            result = lab.validate_claim(report, stage="vs_rho", ledger=self.ledger)
            self.assertEqual(result["status"], "FAIL", f"variant {index}: {result}")
            self.assertTrue(result["validation_errors"] or result["pairing_errors"])

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

    def test_factor_base_complete_passes(self) -> None:
        report = {
            "n_or_bits": 53,
            "factor_base_size_F": 19928,
            "orbit_count_K": 188,
            "dimension_l_or_dim": 188,
            "construction_method": "point_defined_signed_quotient_eta_1_16",
            "materialized": True,
            "construction_wall_ms": 9051.0,
            "retained_bytes": 89887680,
            "fixture_hash": "a",
            "executable_or_source_hash": "b",
            "host_id": "h",
            "resource_caps": {"common_cap_bytes": 17179869184},
            "seeds": {"direct": 1},
            "claim_boundary_non_claims": ["not key recovery"],
        }
        result = lab.validate_claim(report, stage="factor_base", ledger=self.ledger)
        self.assertEqual(result["status"], "PASS", result)


class HelperTests(unittest.TestCase):
    def test_nested_collection_is_rejected_and_exclusive_precompute_is_charged_once(self) -> None:
        row = {
            "fixture_setup_ms": 0.5,
            "collection_ms": 10.0,
            "linear_solve_ms": 2.0,
            "solution_validation_ms": 1.0,
            "reference_validation_ms": 4.0,
            "target_online_wall_ms": 14.5,
            "timing_breakdown_ms": {
                "target_generation": 1.0,
                "packed_verification": 0.5,
            },
            "target_online_phase_ms": {
                "target_query": 1.5,
                "target_pdp": 5.5,
                "target_relation_check": 4.5,
                "target_descent": 0.0,
                "target_recovery_check": 3.0,
            },
        }
        self.assertTrue(lab.exclusive_online_timing_ok(row))
        historical = dict(row)
        historical["target_online_phase_ms"] = dict(row["target_online_phase_ms"], target_pdp=10.0)
        historical["target_online_wall_ms"] = 17.5
        self.assertFalse(lab.exclusive_online_timing_ok(historical))
        self.assertEqual(
            lab.exclusive_precomputation_ms(dict(row, setup_ms=20.0)),
            34.5,
        )
        self.assertIsNone(lab.exclusive_precomputation_ms({"setup_ms": 20.0}))
    def test_independent_replay_pointer_targets_validation_receipt(self) -> None:
        run = lab.REPO / "research/example/autolab/runs/test-run"
        self.assertEqual(
            lab.independent_replay_pointer(run),
            "research/example/autolab/runs/test-run/validation/independent_replay.json",
        )

    def test_single_target_claim_points_to_separate_validation_receipt(self) -> None:
        target = [17, 23]
        scalar = 3
        direct_rows = [
            {
                "kind": "relation_rank_summary",
                "online_target_count": 1,
                "published_q": target,
                "published_fixture_scalar": scalar,
                "status": "SHARED_FACTOR_LOG_ONE_RELATION",
                "uses_retained_factor_logs": True,
                "linear_solution_verified": True,
                "recovered_fixture_scalar": scalar,
                "target_online_phase_ms": {"target_pdp": 1.0},
                "target_online_wall_ms": 1.0,
                "target_trials": 1,
            }
        ]
        rho_rows = [
            {
                "kind": "rho_public_fixture",
                "published_q": target,
                "published_fixture_scalar": scalar,
                "verified": True,
                "walk_ms": 2.0,
                "validation_ms": 0.1,
                "walk_steps": 12,
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "run"
            (run / "artifacts").mkdir(parents=True)
            direct_bin = run / "direct"
            rho_bin = run / "rho"
            direct_bin.write_text("direct")
            rho_bin.write_text("rho")
            claim = lab.draft_vs_rho_claim(
                beat={"timing_class_goal": "single_target_online", "n": 41, "regime": "koblitz"},
                beat_id="koblitz.vs_rho.n41_charged",
                run=run,
                direct_obs={
                    "stdout": "\n".join(json.dumps(row) for row in direct_rows),
                    "exit_code": 0,
                    "whole_process_wall_ms": 4.0,
                },
                rho_obs={
                    "stdout": "\n".join(json.dumps(row) for row in rho_rows),
                    "exit_code": 0,
                    "whole_process_wall_ms": 4.0,
                },
                binaries={"direct": str(direct_bin), "rho": str(rho_bin)},
            )
            self.assertEqual(claim["target_count"], 1)
            accounting = claim["operation_accounting"]
            self.assertEqual(accounting["ic_online_operations"], 1)
            self.assertEqual(accounting["rho_online_operations"], 12)
            self.assertEqual(accounting["rho_per_ic_native_counter"], 12.0)
            self.assertIsNone(accounting["ops_speedup_online"])
            self.assertEqual(lab.operation_accounting_errors(accounting), [])
            self.assertEqual(
                claim["independent_replay_pointer"],
                str(run / "validation/independent_replay.json"),
            )
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
        with self.assertRaisesRegex(lab.AutolabError, 'exactly one online target'):
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

    def test_factor_base_producer_commands_omit_rho(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["koblitz.factor_base.n53"]
        commands = lab.producer_commands(beat, "koblitz.factor_base.n53", 1, binaries=None)
        self.assertIsNone(commands["rho"])
        self.assertIsNone(commands["rho_argv"])
        self.assertIn("koblitz_rank_fixture", commands["direct"])
        self.assertEqual(commands["direct_argv"][1], "53")

    def test_single_target_command_separates_precompute_from_online_target(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["koblitz.vs_rho.n41_charged"]
        commands = lab.producer_commands(beat, "koblitz.vs_rho.n41_charged", 1, binaries=None)
        self.assertEqual(commands["online_target_count"], 1)
        self.assertEqual(commands["precompute_fixtures"], 1)
        self.assertEqual(commands["direct_producer_fixture_count"], 2)
        self.assertEqual(commands["direct_argv"][-1], "2")
        self.assertEqual(commands["rho_argv"][-3], "1")

    def test_vs_rho_batch_launch_is_rejected_before_runner_lock(self) -> None:
        args = type("Args", (), {
            "beat": "koblitz.vs_rho.n41_charged",
            "fixtures": 32,
        })()
        with mock.patch.object(lab, "RunnerLock") as runner_lock:
            with self.assertRaisesRegex(lab.AutolabError, "exactly one online target"):
                lab.launch(args)
            runner_lock.assert_not_called()

    def test_draft_factor_base_claim_from_stdout(self) -> None:
        protocol = lab.load_protocol()
        beat = protocol["beats"]["koblitz.factor_base.n53"]
        stdout = json.dumps(
            {
                "kind": "point_defined_factor_base",
                "evidence_class": "measured_factor_base_construction",
                "n": 53,
                "factor_base_points": 19928,
                "orbit_columns": 188,
                "total_setup_ms": 9000.0,
                "support_payload_lower_bound_bytes": 89887680,
                "support_table_allocated_bytes": 100663296,
                "frobenius_closed": True,
                "negation_closed": True,
                "subgroup_membership_verified": True,
                "pair_index_mode": "signed_quotient",
                "base_hash": "abc",
                "eta": {"numerator": 1, "denominator": 16},
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            (run / "artifacts").mkdir()
            fake_bin = run / "fake_direct"
            fake_bin.write_text("x")
            claim = lab.draft_factor_base_claim(
                beat=beat,
                beat_id="koblitz.factor_base.n53",
                run=run,
                direct_obs={
                    "stdout": stdout,
                    "stderr": "",
                    "exit_code": 0,
                    "whole_process_wall_ms": 12000.0,
                    "children_peak_rss_bytes": 666206208,
                    "seed": 42,
                },
                binaries={"direct": str(fake_bin)},
            )
            self.assertEqual(claim["stage"], "factor_base")
            self.assertEqual(claim["factor_base_size_F"], 19928)
            self.assertEqual(claim["orbit_count_K"], 188)
            self.assertEqual(claim["retained_bytes"], 89887680)
            self.assertEqual(
                claim["independent_replay_pointer"],
                str(run / "validation/independent_replay.json"),
            )
            self.assertTrue(claim["resource_caps"]["under_cap"])
            ledger = lab.load_ledger(protocol)
            result = lab.validate_claim(claim, stage="factor_base", ledger=ledger)
            self.assertEqual(result["status"], "PASS", result)

    def test_claim_check_cli_roundtrip(self) -> None:
        report = {
            "stage": "vs_rho",
            "n": 13,
            "timing_class": "single_target_online",
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
            self.assertEqual(json.loads(out.read_text())["status"], "FAIL")



class OperationAccountingTests(unittest.TestCase):
    def test_every_rung_accounts_and_matches_the_committed_block(self) -> None:
        for name, cfg in op_accounting.RUNGS.items():
            block = op_accounting.account_rung(name)
            self.assertEqual(lab.operation_accounting_errors(block), [], name)
            claim = json.loads((lab.REPO / cfg["dir"] / "claim_report_vs_rho.json").read_text())
            self.assertEqual(claim.get("operation_accounting"), block, f"{name}: rerun op_accounting.py --write")
            self.assertIsNone(block["total_work_S"], name)
            self.assertIsNone(block["ops_speedup_online"], name)

    def test_n73_separates_contended_runs_and_native_counter_ratios(self) -> None:
        block = op_accounting.account_rung("n73")
        self.assertEqual(block["ic_online_operations"], 457561)
        self.assertEqual(block["host_contention"]["contended_runs"], ["R2", "R3"])
        self.assertAlmostEqual(block["wall_speedup_median_uncontended_runs"], 184.8588, places=3)
        self.assertGreater(block["rank_mean_to_frozen_probe_ratio"], 100)
        self.assertIn("not an estimate of unseen-target cost", block["ic_rank_probes_note"])
        self.assertLess(block["rho_expected_steps_per_rank_mean_probe"], 1)
        expected = op_accounting.expected_rho_steps(86020738150056119, 73)
        self.assertAlmostEqual(block["rho_per_ic_native_counter"], expected / 457561)

    def test_native_counters_do_not_promote_total_work(self) -> None:
        block = op_accounting.account_rung("n73")
        self.assertEqual(block["comparison_status"], "native_counters_only")
        self.assertIsNone(block["total_work_S"])
        self.assertIsNone(block["ops_speedup_online"])


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


class BatchPanelTests(unittest.TestCase):
    BEAT = "koblitz.compact_orbit.n61_panel"

    def setUp(self) -> None:
        self.protocol = lab.load_protocol()
        self.ledger = lab.load_ledger(self.protocol)
        self.beat = self.protocol["beats"][self.BEAT]

    def test_panel_beat_is_registered_as_diagnostic_only(self) -> None:
        self.assertEqual(self.beat["launch_mode"], "batch_panel")
        self.assertEqual(self.beat["stage"], "end_to_end_dlp")
        self.assertFalse(self.beat["primary_speedup_eligible"])
        self.assertGreaterEqual(self.beat["targets_minimum"], 1024)
        for targets in self.beat["k_candidates"]:
            self.assertGreaterEqual(int(targets), self.beat["targets_minimum"])
        rows = {row["beat_id"]: row for row in lab.plan(self.protocol, self.ledger)["beats"]}
        self.assertFalse(rows[self.BEAT]["primary_speedup_eligible"])

    def test_single_target_launch_rejects_the_panel_beat(self) -> None:
        args = lab.parser().parse_args(["launch", "--beat", self.BEAT])
        with self.assertRaisesRegex(lab.AutolabError, "launch-panel"):
            lab.launch(args)

    def test_panel_rejects_small_batches_and_single_target_beats(self) -> None:
        args = lab.parser().parse_args(["launch-panel", "--beat", self.BEAT, "--targets", "32"])
        with self.assertRaisesRegex(lab.AutolabError, ">= 1024"):
            lab.launch_panel(args)
        args = lab.parser().parse_args(["launch-panel", "--beat", "koblitz.vs_rho.n37_wall"])
        with self.assertRaisesRegex(lab.AutolabError, "not a batch panel"):
            lab.launch_panel(args)

    def test_parse_time_l_ignores_producer_lines(self) -> None:
        text = textwrap.dedent(
            """\
            note: producer line 12
                  295.04 real       273.93 user        16.72 sys
               11502862336  maximum resident set size
                  2274502049018  instructions retired
                       12  swaps
            """
        )
        parsed = lab.parse_time_l(text)
        self.assertEqual(parsed["max_rss_bytes"], 11502862336)
        self.assertEqual(parsed["instructions_retired"], 2274502049018)
        self.assertEqual(parsed["swaps"], 12)
        self.assertNotIn("page_faults", parsed)

    def block_records(self, scalars):
        ic = [
            {"kind": "compact_orbit_dlp_target", "fixture_index": i, "published_fixture_scalar": s,
             "target": [s, s + 1], "recovered_matches_published": True, "group_verified": True,
             "probes": 3, "target_ms": 0.5}
            for i, s in enumerate(scalars)
        ]
        rho = [
            {"kind": "rho_ks_batch_fixture", "fixture_index": i, "published_fixture_scalar": s,
             "recovered_fixture_scalar": s, "published_q": [s, s + 1], "verified": True}
            for i, s in enumerate(scalars)
        ] + [{"kind": "rho_ks_batch_summary"}]
        return ic, rho

    def test_block_check_requires_same_verified_targets(self) -> None:
        ic, rho = self.block_records([5, 9])
        checks = lab.check_panel_block(ic, rho, [5, 9])
        for key in ("ic_all_verified", "rho_all_verified", "ic_matches_corpus",
                    "rho_matches_corpus", "same_target_points"):
            self.assertTrue(checks[key], key)
        rho[1]["published_q"] = [9, 11]
        self.assertFalse(lab.check_panel_block(ic, rho, [5, 9])["same_target_points"])
        ic[0]["group_verified"] = False
        self.assertFalse(lab.check_panel_block(ic, rho, [5, 9])["ic_all_verified"])
        self.assertFalse(lab.check_panel_block([], rho, [5, 9])["ic_all_verified"])

    def test_untimed_digest_ignores_timers(self) -> None:
        ic, _ = self.block_records([5])
        slower = [dict(ic[0], target_ms=99.0)]
        self.assertEqual(lab.untimed_digest(ic), lab.untimed_digest(slower))

    def test_panel_claims_pass_end_to_end_and_fail_vs_rho_closed(self) -> None:
        with tempfile.TemporaryDirectory(dir=lab.RUNS_DIR.parent) as tmp:
            run = Path(tmp)
            (run / "artifacts").mkdir()
            (run / "artifacts/replay_ic.json").write_text("{}\n")
            summary = {
                "targets": 1024, "K": 600, "subgroup_order_bits": 47.21,
                "corpora": {"eval": {"name": "e"}, "tune": {"name": "t"}, "disjoint": True},
                "executables": {"git_head": "0" * 40},
                "host": {"node": "test"},
                "isolation": "none",
                "wall_ratio": {"median": 0.3}, "user_ratio": None, "instructions_ratio": None,
                "verification": {"ic_all_verified": True, "rho_all_verified": True},
                "blocks": [{"block": 0, "ic_timing_ms": {"process_total": 1.0}}],
            }
            end_to_end, vs_rho = lab.draft_panel_claims(
                beat_id=self.BEAT, beat=self.beat, run=run, summary=summary
            )
        e2e = lab.validate_claim(end_to_end, stage="end_to_end_dlp", ledger=self.ledger)
        self.assertEqual(e2e["status"], "PASS", e2e)
        check = lab.validate_claim(vs_rho, stage="vs_rho", ledger=self.ledger)
        self.assertEqual(check["status"], "FAIL")
        self.assertIn("target_count must equal 1", check["validation_errors"])
        for absent in ("ic_online_wall_ms", "rho_online_wall_ms", "ic_online_phase_ms", "online_interval"):
            self.assertNotIn(absent, vs_rho)

    def test_rss_model_covers_measured_peaks_and_the_table_doubling(self) -> None:
        # Measured n=61 IC peaks (GiB) from the 1,024/4,096-target tunes, 2026-10-01.
        measured = {400: 1.23, 500: 1.35, 600: 2.51, 700: 2.68, 800: 4.89, 1000: 5.39, 1200: 10.00}
        for k, gib in measured.items():
            ratio = lab.estimate_ic_rss_bytes(self.beat, k) / 2**30 / gib
            self.assertTrue(1.0 <= ratio <= 1.05, (k, ratio))
        self.assertLess(lab.estimate_ic_rss_bytes(self.beat, 1480), 12 * 2**30)
        self.assertGreater(lab.estimate_ic_rss_bytes(self.beat, 1500), 19 * 2**30)
        legacy = {key: value for key, value in self.beat.items() if key != "ic_rss_model"}
        legacy["ic_bytes_per_regular_state"] = 94
        self.assertEqual(lab.estimate_ic_rss_bytes(legacy, 1000), 94 * 61 * 1000 * 1000)

    def test_tune_targets_option_is_bounded_and_has_grids(self) -> None:
        args = lab.parser().parse_args(
            ["launch-panel", "--beat", self.BEAT, "--targets", "4096", "--tune-targets", "512"]
        )
        self.assertEqual(args.tune_targets, 512)
        self.assertFalse(args.tune_only)
        with self.assertRaisesRegex(lab.AutolabError, "tune targets must be >= 1024"):
            lab.launch_panel(args)
        self.assertIsNone(lab.parser().parse_args(["launch-panel", "--beat", self.BEAT]).tune_targets)
        grids = self.beat["k_candidates_by_tune_targets"]
        for tune_targets, by_l in grids.items():
            self.assertGreaterEqual(int(tune_targets), self.beat["targets_minimum"])
            for targets, ks in by_l.items():
                self.assertIn(targets, self.beat["k_candidates"])
                self.assertNotEqual(targets, tune_targets)
                for k in ks:
                    self.assertLess(lab.estimate_ic_rss_bytes(self.beat, k), 12 * 2**30)
        tune = lab.panel_corpus_name(self.beat["corpora"]["tune"], 1024)
        for targets in grids["1024"]:
            self.assertNotEqual(tune, lab.panel_corpus_name(self.beat["corpora"]["eval"], int(targets)))

    def test_panel_manifest_never_lists_itself(self) -> None:
        with tempfile.TemporaryDirectory(dir=lab.RUNS_DIR) as tmp:
            run = Path(tmp)
            (run / "artifacts").mkdir()
            (run / "artifacts/k_tune.json").write_text("{}\n")
            lab.write_panel_manifest(run)
            (run / "artifacts/blocks.json").write_text("[]\n")
            lab.write_panel_manifest(run)
            files = lab.read_json(run / "artifacts/review_manifest.json")["files"]
            self.assertEqual(sorted(files), ["artifacts/blocks.json", "artifacts/k_tune.json"])
            result = lab.verify(lab.parser().parse_args(["verify", "--run-id", run.name]))
            self.assertEqual(result["status"], "PASS", result)

    def test_panel_resume_needs_an_existing_run(self) -> None:
        args = lab.parser().parse_args(["launch-panel", "--beat", self.BEAT, "--resume", "no-such-run"])
        with self.assertRaisesRegex(lab.AutolabError, "no run to resume"):
            lab.launch_panel(args)


class SingleTargetPanelTests(unittest.TestCase):
    BEAT = "koblitz.compact_orbit.n61_single_target"

    @classmethod
    def setUpClass(cls) -> None:
        import single_target_panel

        cls.stp = single_target_panel
        cls.protocol = lab.load_protocol()
        cls.ledger = lab.load_ledger(cls.protocol)
        cls.beat = cls.protocol["beats"][cls.BEAT]
        cls.curve = single_target_panel.cached_curve(cls.beat["curve_fixture"])

    def test_beat_is_one_target_against_strong_rho(self) -> None:
        self.assertEqual(self.beat["launch_mode"], "single_target_panel")
        self.assertEqual(self.beat["stage"], "vs_rho")
        self.assertEqual(self.beat["panel_producers"]["rho"]["env"]["KIC_RHO_RUNG"], "3")
        for producer in list(self.beat["panel_producers"].values()) + list(self.beat["frozen_originals"].values()):
            self.assertTrue((lab.REPO / producer["source"]).is_file(), producer["source"])
        self.assertNotEqual(self.beat["panel_producers"]["ic"]["source"], self.beat["frozen_originals"]["ic"]["source"])
        rows = {row["beat_id"]: row for row in lab.plan(self.protocol, self.ledger)["beats"]}
        self.assertFalse(rows[self.BEAT]["primary_speedup_eligible"])

    def test_other_launchers_reject_the_beat(self) -> None:
        with self.assertRaisesRegex(lab.AutolabError, "launch-single"):
            lab.launch(lab.parser().parse_args(["launch", "--beat", self.BEAT]))
        with self.assertRaisesRegex(lab.AutolabError, "not a batch panel"):
            lab.launch_panel(lab.parser().parse_args(["launch-panel", "--beat", self.BEAT]))
        args = lab.parser().parse_args(["launch-single", "--beat", self.BEAT, "--workloads", "4"])
        with self.assertRaisesRegex(lab.AutolabError, ">= 16"):
            self.stp.launch_single(args, lab)

    def test_public_points_are_deterministic_subgroup_points(self) -> None:
        domain = self.beat["target_law"]["domain"]
        first = self.stp.public_point(self.curve, domain, "eval", 0)
        self.assertEqual(first, self.stp.public_point(self.curve, domain, "eval", 0))
        self.assertNotEqual(first, self.stp.public_point(self.curve, domain, "tune", 0))
        self.assertIsNone(self.curve.mul(first, self.curve.r))
        self.assertEqual(self.curve.decode(list(first)), first)

    def test_replay_checks_digest_scalar_and_relation(self) -> None:
        scalar = 987654321
        target = list(self.curve.mul(self.curve.g, scalar))
        fixture = dict(self.beat["curve_fixture"], targets=[[str(v) for v in target]])
        generator = [int(v) for v in fixture["generator"]]
        cert = {"arm": "rho", "curve_id": "c", "generator": generator, "target": target, "scalar": scalar}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cert.json"
            path.write_bytes(self.stp.identity.canonical(cert) + b"\n")
            digest = self.stp.identity.sha256(cert)
            self.assertTrue(self.stp.replay_certificate(path, digest, fixture, None)["statement_holds"])
            self.assertFalse(self.stp.replay_certificate(path, "0" * 64, fixture, None)["statement_holds"])
            points = [self.curve.mul(self.curve.g, s) for s in (11, 22, 33)]
            last = self.curve.add(self.curve.neg(self.curve.add(self.curve.add(points[0], points[1]), points[2])),
                                  tuple(target))
            base = points + [last]
            ic = {"arm": "ic", "curve_id": "c", "generator": generator, "target": target, "scalar": scalar,
                  "relation": {"base_hash": "b", "point_indices": [0, 1, 2, 3], "x_codes": [p[0] for p in base]}}
            path.write_bytes(self.stp.identity.canonical(ic) + b"\n")
            replay = self.stp.replay_certificate(path, self.stp.identity.sha256(ic), fixture, base)
            self.assertTrue(replay["statement_holds"], replay)
            replay = self.stp.replay_certificate(path, self.stp.identity.sha256(ic), fixture, base[:3] + [points[0]])
            self.assertFalse(replay["checks"]["relation_sums_to_target"])

    def test_untimed_view_drops_timers_only(self) -> None:
        record = {"kind": "k", "probes": 3, "online_ms": 1.0, "target_generation_ms_excluded": 0.1,
                  "online_start_ns": 5, "online_stop_event": "e", "relation_checks": 1}
        self.assertEqual(self.stp.untimed(record, ("relation_checks",)), {"kind": "k", "probes": 3})


if __name__ == "__main__":
    unittest.main()
