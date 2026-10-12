"""Source-only controls for the proposed retained-N83 mixed-orbit gate."""

import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import mixed_orbit_supervisor as gate


def config() -> dict:
    return {
        "schema": "n83.mixed-orbit-partial-config/v1", "study": gate.STUDY,
        "small_object_blake3": gate.SMALL_BLAKE3,
        "envelope_object_blake3": gate.ENVELOPE_BLAKE3,
        "public_corpus_path": "/corpus.json",
        "public_corpus_canonical_json_blake3": gate.CORPUS_BLAKE3,
        "budget_seconds": gate.WORKER_WALL, "max_trials": gate.TRIALS,
        "conflict_budget_per_trial": gate.CONFLICTS,
        "max_models_per_trial": gate.MODELS,
        "small_summands": gate.SMALL_SUMMANDS,
        "residual_summands": gate.RESIDUAL_SUMMANDS,
        "order_first_pair": True, "max_variables": 150000,
        "max_domain_clauses": 12000000,
        "memory_cgroup_limit_bytes": 4 * 1024**3,
        "memory_cgroup_swap_limit_bytes": 0,
        "graph_limits": {"max_merge_depth": 32, "max_row_terms": 4096,
                         "max_pivots": 1000000},
        "source_blake3": "a" * 64,
    }


def result() -> dict:
    return {
        "schema": "n83.mixed-orbit-partial-result/v1", "study": gate.STUDY,
        "status": "UNKNOWN_search_cap", "small_object_blake3": gate.SMALL_BLAKE3,
        "envelope_object_blake3": gate.ENVELOPE_BLAKE3,
        "public_corpus_canonical_json_blake3": gate.CORPUS_BLAKE3,
        "small_signed_orbits": 1182, "residual_signed_orbits": 866,
        "trials_attempted": 1, "solver_exhausted_trials": 1,
        "solver_refuted_trials": 0, "completed_rows": 0,
        "solver_calls": 1, "solver_conflicts": 5000,
        "solver_models": 0, "invalid_models": 0,
        "group_rejected_models": 0, "graph_accepted_candidates": 0,
        "graph_pivot_rows": 0, "graph_full_from_partials": 0,
        "graph_capped_rows": 0, "graph_pivot_reductions": 0,
        "base_import_elapsed_ms_L0": 9.0, "worker_wall_ms_L0": 30.0,
        "memory_cgroup_peak_bytes": 250_000_000,
        "matrix_rank": None, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


class MixedOrbitSupervisorTests(unittest.TestCase):
    def test_additional_budget_is_required_before_build(self):
        with tempfile.TemporaryDirectory(prefix="n83-mixed-budget-") as folder:
            root = Path(folder)
            build, output = root / "build", root / "run"
            auth = root / "authorization.json"
            record = {
                "schema": "n83.additional-local-budget/v2", "study": gate.STUDY,
                "purpose": "mixed_orbit_two_residual_first_gate",
                "prior_total_charged_seconds": "12600",
                "prior_import_audit_sha256": gate.PRIOR_IMPORT_AUDIT_SHA256,
                "build_dir": str(build), "output_dir": str(output),
                "additional_authorized_seconds": "2399",
                "user_approval_reference": "synthetic-test-only",
            }
            auth.write_text(json.dumps(record))
            with patch.object(gate, "image_identity",
                              side_effect=AssertionError("Docker must not be inspected")):
                with self.assertRaisesRegex(ValueError, "below full conservative charge"):
                    gate.build(auth, build, output)
            self.assertFalse(build.exists())
            record["additional_authorized_seconds"] = "2400"
            self.assertIsNone(gate.authorization_error(record, build, output))
            record["prior_total_charged_seconds"] = "11700"
            self.assertIn("prior_total", gate.authorization_error(record, build, output))
            record["prior_total_charged_seconds"] = "12600"
            record["prior_import_audit_sha256"] = "0" * 64
            self.assertIn("prior_import_audit", gate.authorization_error(record, build, output))
            record["prior_import_audit_sha256"] = gate.PRIOR_IMPORT_AUDIT_SHA256
            record["schema"] = "n83.additional-local-budget/v1"
            self.assertIn("schema", gate.authorization_error(record, build, output))

    def test_prior_import_receipt_is_exact_and_charged(self):
        self.assertIsNone(gate.prior_import_evidence_error(gate.ROOT))
        with tempfile.TemporaryDirectory(prefix="n83-mixed-prior-") as folder:
            root = Path(folder)
            path = root / gate.PRIOR_IMPORT_AUDIT
            path.parent.mkdir(parents=True)
            path.write_text(json.dumps({"status": "PASS", "conservative_charge_seconds": 900}))
            self.assertIn("bytes changed", gate.prior_import_evidence_error(root))
            with patch.object(gate, "sha256_file",
                              return_value=gate.PRIOR_IMPORT_AUDIT_SHA256):
                self.assertIn("disposition", gate.prior_import_evidence_error(root))

    def test_receipts_reject_false_success_and_bad_counters(self):
        verify = {
            "schema": "n83.mixed-orbit-import/v1", "status": "PASS",
            "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
            "small_object_blake3": gate.SMALL_BLAKE3,
            "envelope_object_blake3": gate.ENVELOPE_BLAKE3,
            "small_points": 196212, "small_signed_orbits": 1182,
            "residual_points": 143756, "residual_signed_orbits": 866,
            "solver_search_executed": False,
        }
        self.assertIsNone(gate.verify_receipt_error(verify))
        verify["residual_signed_orbits"] -= 1
        self.assertIn("residual_signed_orbits", gate.verify_receipt_error(verify))
        c, r = config(), result()
        self.assertIsNone(gate.search_receipt_error(c, r))
        r["total_index_calculus_runtime_ms"] = 30.0
        self.assertIn("total_index_calculus_runtime_ms", gate.search_receipt_error(c, r))
        r = result()
        r["status"] = "PASS_bounded_partial_search"
        self.assertIn("hides a cap", gate.search_receipt_error(c, r))
        r = result()
        r["graph_accepted_candidates"] = 1
        self.assertIn("counters disagree", gate.search_receipt_error(c, r))
        cap = {"schema": "n83.mixed-orbit-partial-cap/v1",
               "status": "UNKNOWN_budget", "budget_seconds": gate.WORKER_WALL,
               "elapsed_process_wall_ms": 1150000.0}
        self.assertTrue(gate.valid_worker_cap(cap))
        cap["budget_seconds"] -= 1
        self.assertFalse(gate.valid_worker_cap(cap))

    def test_build_and_run_preserve_pinned_stage_commands(self):
        with tempfile.TemporaryDirectory(prefix="n83-mixed-stages-") as folder:
            root = Path(folder)
            build_dir, output = root / "build", root / "run"
            cargo_home = root / "cargo"
            (cargo_home / "registry").mkdir(parents=True)
            small, envelope = root / "small.gz", root / "envelope.gz"
            small.write_bytes(b"synthetic small object")
            envelope.write_bytes(b"synthetic envelope object")
            auth = root / "authorization.json"
            auth.write_text(json.dumps({
                "schema": "n83.additional-local-budget/v2", "study": gate.STUDY,
                "purpose": "mixed_orbit_two_residual_first_gate",
                "prior_total_charged_seconds": "12600",
                "prior_import_audit_sha256": gate.PRIOR_IMPORT_AUDIT_SHA256,
                "build_dir": str(build_dir), "output_dir": str(output),
                "additional_authorized_seconds": "2400",
                "user_approval_reference": "synthetic-test-only",
            }))
            real_sha256 = gate.sha256_file

            def object_hash(path):
                if Path(path) == small:
                    return gate.SMALL_SHA256
                if Path(path) == envelope:
                    return gate.ENVELOPE_SHA256
                return real_sha256(path)

            launched = []

            def fake_stage(*, name, command, wall, output, container_name):
                launched.append((name, command, wall))
                self.assertEqual(command[command.index("--network") + 1], "none")
                self.assertEqual(command[command.index("--pull") + 1], "never")
                if name == "build":
                    self.assertEqual(wall, gate.BUILD_WALL)
                    artifact = output / "target/release/examples/koblitz_n83_mixed_orbit_worker"
                    artifact.parent.mkdir(parents=True)
                    artifact.write_bytes(b"\x7fELFsynthetic binary")
                else:
                    self.assertEqual(command[command.index("--memory") + 1], "4096m")
                    self.assertEqual(command[command.index("--memory-swap") + 1], "4096m")
                    self.assertEqual(command[command.index("--cpus") + 1], "1")
                    if name == "verify":
                        self.assertEqual(wall, gate.VERIFY_CAP)
                        (output / "verify.stdout").write_text(json.dumps({
                            "schema": "n83.mixed-orbit-import/v1", "status": "PASS",
                            "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
                            "small_object_blake3": gate.SMALL_BLAKE3,
                            "envelope_object_blake3": gate.ENVELOPE_BLAKE3,
                            "small_points": 196212, "small_signed_orbits": 1182,
                            "residual_points": 143756, "residual_signed_orbits": 866,
                            "solver_search_executed": False,
                        }))
                    else:
                        self.assertEqual(wall, gate.SEARCH_CAP)
                        (output / "run").mkdir()
                        (output / "run/config.json").write_text(json.dumps(config()))
                        (output / "run/result.json").write_text(json.dumps(result()))
                        (output / "run/events.jsonl").write_text(json.dumps({
                            "trial": 0, "a": gate.FIRST_TRIAL_A,
                            "b": gate.FIRST_TRIAL_B,
                            "elapsed_process_wall_ms": 10.0,
                            "search": {"solver_calls": 1, "conflicts": 5000,
                                       "models": 0, "invalid_models": 0,
                                       "group_rejected_models": 0,
                                       "exhausted": True, "refuted": False},
                            "candidate": None, "graph": None,
                        }) + "\n")
                return {"name": name, "exit_code": 0, "timed_out": False,
                        "cleanup_verified": True, "launch_error": None,
                        "container_state": {"ExitCode": 0, "OOMKilled": False}}

            with patch.dict(os.environ, {"CARGO_HOME": str(cargo_home)}), \
                 patch.object(gate, "clean_commit", return_value="a" * 40), \
                 patch.object(gate, "sha256_file", side_effect=object_hash), \
                 patch.object(gate, "image_identity", return_value=gate.IMAGE), \
                 patch.object(gate, "stage", side_effect=fake_stage):
                built = gate.build(auth, build_dir, output)
                self.assertEqual(built["status"], "PASS_build_only")
                outer = gate.run(auth, build_dir, output, small, envelope)
            self.assertEqual(outer["status"], "UNKNOWN_search_cap")
            self.assertEqual([name for name, _, _ in launched], ["build", "verify", "search"])
            self.assertEqual(json.loads((output / "budget.json").read_text())
                             ["conservative_charge_seconds"], 2400)

    def test_independent_partial_and_completed_row_replay(self):
        generator = gate.PRIMARY_GENERATOR
        negative = (generator[0], generator[0] ^ generator[1])
        self.assertTrue(gate.replay_equation("1", "0", [generator] * 3,
                                             [negative] * 2, generator))
        self.assertFalse(gate.replay_equation("1", "0", [generator] * 3,
                                              [negative, generator], generator))
        a, b = int(gate.FIRST_TRIAL_A), int(gate.FIRST_TRIAL_B)
        target_scalar = ((1 - a) * pow(b, -1, gate.PRIMARY_ORDER)) % gate.PRIMARY_ORDER
        target = gate.point_mul(generator, target_scalar)
        self.assertIsNotNone(target)
        with tempfile.TemporaryDirectory(prefix="n83-mixed-replay-") as folder:
            root = Path(folder)
            archive = root / "small.jsonl.gz"
            envelope = root / "envelope.jsonl.gz"
            point = [str(generator[0]), str(generator[1])]
            negative_point = [str(negative[0]), str(negative[1])]
            twice = gate.point_mul(generator, 2)
            negative_twice = (twice[0], twice[0] ^ twice[1])
            negative_twice_point = [str(negative_twice[0]), str(negative_twice[1])]
            with gzip.open(archive, "wt", encoding="utf-8") as stream:
                stream.write(json.dumps({"representatives": [point] * 1182,
                                         "point_count": 3}) + "\n")
                for _ in range(3):
                    stream.write(json.dumps({"point": point}) + "\n")
            with gzip.open(envelope, "wt", encoding="utf-8") as stream:
                stream.write(json.dumps({"point_count": 5}) + "\n")
                for _ in range(3):
                    stream.write(json.dumps({"point": point}) + "\n")
                stream.write(json.dumps({"point": negative_point}) + "\n")
                stream.write(json.dumps({"point": negative_twice_point}) + "\n")
            event = {
                "trial": 0, "a": gate.FIRST_TRIAL_A,
                "b": gate.FIRST_TRIAL_B, "elapsed_process_wall_ms": 1.0,
                "search": {"solver_calls": 1, "conflicts": 20, "models": 1,
                           "invalid_models": 0, "group_rejected_models": 0,
                           "exhausted": False, "refuted": False},
                "candidate": {
                    "source_id": 1, "a": gate.FIRST_TRIAL_A,
                    "b": gate.FIRST_TRIAL_B,
                    "small_point_indices": [0, 1, 2],
                    "residual_points": [negative_point] * 2,
                },
                "graph": {"kind": "full", "a": gate.FIRST_TRIAL_A,
                          "b": gate.FIRST_TRIAL_B,
                          "small": [["0", "1"]], "source_coeffs": [["1", "1"]],
                          "merge_depth": 0},
            }
            events = root / "events.jsonl"
            events.write_text(json.dumps(event) + "\n")
            report = result()
            report.update({
                "status": "PASS_bounded_partial_search", "solver_exhausted_trials": 0,
                "solver_conflicts": 20, "solver_models": 1,
                "graph_accepted_candidates": 1, "completed_rows": 1,
            })
            with patch.object(gate, "SMALL_POINTS", 3), \
                 patch.object(gate, "ENVELOPE_POINTS", 5):
                self.assertIsNone(gate.replay_events(events, archive, envelope, target, report))
                event["a"] = "1"
                events.write_text(json.dumps(event) + "\n")
                self.assertIn("frozen first-trial", gate.replay_events(
                    events, archive, envelope, target, report))
                event["a"] = gate.FIRST_TRIAL_A
                event["candidate"]["residual_points"][0] = negative_twice_point
                events.write_text(json.dumps(event) + "\n")
                self.assertIn("partial point equation", gate.replay_events(
                    events, archive, envelope, target, report))
                triple = gate.point_mul(generator, 3)
                negative_thrice = (triple[0], triple[0] ^ triple[1])
                outside = [str(negative_thrice[0]), str(negative_thrice[1])]
                self.assertTrue(gate.replay_equation(
                    gate.FIRST_TRIAL_A, gate.FIRST_TRIAL_B,
                    [generator] * 3, [generator, negative_thrice], target))
                event["candidate"]["residual_points"] = [point, outside]
                events.write_text(json.dumps(event) + "\n")
                self.assertIn("outside pinned envelope suffix", gate.replay_events(
                    events, archive, envelope, target, report))
                event["candidate"]["residual_points"] = [negative_point] * 2
                event["graph"]["source_coeffs"] = [["1", "2"]]
                events.write_text(json.dumps(event) + "\n")
                self.assertIn("source combination", gate.replay_events(
                    events, archive, envelope, target, report))
                event["graph"]["source_coeffs"] = [["1", "1"]]
                event["graph"]["small"] = [["0", "2"]]
                events.write_text(json.dumps(event) + "\n")
                self.assertIn("completed row failed", gate.replay_events(
                    events, archive, envelope, target, report))


if __name__ == "__main__":
    unittest.main()
