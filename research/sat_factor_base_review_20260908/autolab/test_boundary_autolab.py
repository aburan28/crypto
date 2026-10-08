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
