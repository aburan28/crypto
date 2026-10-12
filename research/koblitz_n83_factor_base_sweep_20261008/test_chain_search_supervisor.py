"""Synthetic process guards for retained N83 chained-S3 search."""

import copy
import gzip
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import chain_search_supervisor as guard


def receipts(config: dict, status: str = "PREFLIGHT_ONLY") -> tuple[dict, dict]:
    limits = {
        "max_variables": config["max_variables"],
        "max_domain_clauses": config["max_domain_clauses"],
        "max_models": config["max_models"],
        "conflict_budget": config["conflict_budget"],
    }
    trials = 0 if status == "PREFLIGHT_ONLY" else 1
    calls = 0 if status == "PREFLIGHT_ONLY" else 1
    verified = status == "PASS_verified_target_only"
    unknown = status == "UNKNOWN_solver_cap"
    report = {key: 0 for key in guard.COUNTERS}
    report.update({
        "factor_base_size": 166 * config["columns"],
        "orbit_count": config["columns"], "trials": trials,
        "sat_calls": calls, "sat_unknowns": int(unknown),
        "sat_models": int(verified),
        "relations": int(verified), "independent_relations": int(verified),
        "linear_solve_attempts": int(verified),
        "matrix_rows": int(verified), "matrix_columns": config["columns"] + 1,
        "terminal_matrix_rank": int(verified), "rank_checks": int(verified),
        "relation_collection_ns": "1000", "linear_algebra_ns": "0",
    })
    worker_config = {
        "schema": guard.WORKER_CONFIG_SCHEMA, "study": guard.STUDY,
        "curve_a": 0, "fixture": 0, "orbit_columns": config["columns"],
        "summands": config["summands"], "strategy": guard.strategy_name(config),
        "policy": config["policy"], "seed": config["seed"],
        "max_trials": config["max_trials"],
        "budget_seconds": config["worker_wall_seconds"],
        "allow_direct_relation": False,
        "source_commit": config["source_commit"],
        "source_attestation_mode": "supervised_snapshot",
        "fixture_dir": "/fixtures",
        "chain_s3_limits": limits,
        "memory_cgroup_limit_bytes": config["memory_cgroup_limit_bytes"],
        "memory_cgroup_swap_limit_bytes": 0,
        "source_adapter_blake3": "a" * 64,
    }
    summary = {
        **worker_config, "schema": guard.WORKER_SCHEMA,
        "status": status, "object": config["selected_object"],
        "point_set_blake3": config["selected_point_set_blake3"],
        "panel_manifest_blake3": "d" * 64,
        "public_corpus_canonical_json_blake3": "e" * 64,
        "source_chain_blake3": "f" * 64,
        "source_index_calculus_blake3": "1" * 64,
        "memory_cgroup_peak_bytes": 12345,
        "base_import_ms": 1.0, "target_validation_ms": 2.0,
        "solver_ms": 3.0, "post_solver_validation_ms": 0.0,
        "pre_summary_process_wall_ms": 6.0,
        "column_log_verification": False,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
        "solver_stage_executed": calls > 0,
        "report": report,
        "attempt_records": ([{
            "trial": 1, "a": "1", "b": "2", "target": ["3", "4"],
            "disposition": "RelationFound" if verified else "Unknown",
            "decomposition_indices": [0] * config["summands"] if verified else None,
            "solver_calls": 1, "models": int(verified), "conflicts": 0,
            "implied_rows": int(verified),
        }] if calls else []),
        "relation_matrix": ([{
            "a": "1", "b": "2", "summands": [[0, 0]] * config["summands"],
            "summand_negated": [False] * config["summands"],
            "row": [str(config["summands"])] + ["0"] * (config["columns"] - 1),
        }] if verified else []),
        "rank_history": ([{
            "rows": 1, "columns": config["columns"] + 1, "rank": 1,
            "candidate_produced": True, "candidate_verified": True,
        }] if verified else []),
        "verified_log": "3" if verified else None,
    }
    return worker_config, summary


class ChainSearchGuardTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "columns": 64, "policy": "public_x_hash", "seed": guard.SEEDS[0],
            "summands": 5, "max_trials": 0, "max_variables": 100000,
            "max_domain_clauses": 300000, "max_models": 1,
            "conflict_budget": 1000, "worker_wall_seconds": 1,
            "memory_cgroup_limit_bytes": 256 * 1024 * 1024,
            "source_commit": "a" * 40,
            "selected_object": "objects/" + "b" * 64 + ".jsonl.gz",
            "selected_point_set_blake3": "c" * 64,
        }

    def test_full_augmented_rank_is_admitted(self):
        config = dict(self.config, max_trials=65)
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        # Sixty-four orbit pivots and one target-column pivot are possible.
        # This checks receipt structure; group replay is checked separately.
        attempt = worker["attempt_records"][0]
        relation = worker["relation_matrix"][0]
        attempts, relations = [], []
        for i in range(65):
            column = min(i, 63)
            next_attempt = copy.deepcopy(attempt)
            next_attempt["trial"] = i + 1
            next_attempt["b"] = "2" if i == 64 else "1"
            next_attempt["decomposition_indices"] = [166 * column] * config["summands"]
            attempts.append(next_attempt)
            next_relation = copy.deepcopy(relation)
            next_relation["b"] = next_attempt["b"]
            next_relation["summands"] = [[column, 0]] * config["summands"]
            next_relation["row"] = ["0"] * config["columns"]
            next_relation["row"][column] = str(config["summands"])
            relations.append(next_relation)
        worker["attempt_records"] = attempts
        worker["relation_matrix"] = relations
        worker["rank_history"][0].update(rows=65, rank=65)
        worker["report"].update(
            trials=65, relations=65, independent_relations=65,
            sat_calls=65, sat_models=65, matrix_rows=65,
            terminal_matrix_rank=65,
        )
        self.assertIsNone(guard.receipt_error(worker, worker_config, config))
        worker["report"]["independent_relations"] = 66
        self.assertIn("rank", guard.receipt_error(worker, worker_config, config))
        worker["report"]["independent_relations"] = 64
        self.assertIn("rank", guard.receipt_error(worker, worker_config, config))

    def test_receipt_checks_rank_solver_caps_and_false_completion(self):
        worker_config, worker = receipts(self.config)
        self.assertIsNone(guard.receipt_error(worker, worker_config, self.config))
        worker["report"]["linear_solve_attempts"] = 1
        worker["report"]["rank_checks"] = 1
        worker["rank_history"] = [{
            "rows": 0, "columns": self.config["columns"] + 1, "rank": 0,
            "candidate_produced": False, "candidate_verified": False,
        }]
        self.assertIsNone(guard.receipt_error(worker, worker_config, self.config))
        worker["report"]["relations"] = 1
        worker["report"]["independent_relations"] = 1
        self.assertIsNotNone(guard.receipt_error(worker, worker_config, self.config))
        config = dict(self.config, max_trials=1)
        worker_config, worker = receipts(config, "UNKNOWN_solver_cap")
        self.assertIsNone(guard.receipt_error(worker, worker_config, config))
        worker["report"]["sat_unknowns"] = 0
        self.assertIn("solver-cap", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        self.assertIsNone(guard.receipt_error(worker, worker_config, config))
        worker["relation_matrix"][0]["row"].pop()
        self.assertIn("relation replay", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["relation_matrix"][0]["row"][0] = str(config["summands"] + 1)
        self.assertIn("signed-Frobenius", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["attempt_records"][0]["decomposition_indices"][0] = 1
        self.assertIn("labels", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["relation_matrix"][0]["a"] = "7"
        self.assertIn("labels", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["rank_history"].clear()
        self.assertIn("replay ledger", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["attempt_records"][0]["solver_calls"] = 0
        self.assertIn("SAT counters", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["report"]["independent_relations"] = 0
        self.assertIsNotNone(guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["total_index_calculus_runtime_ms"] = 6.0
        self.assertIn("unmeasured", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker_config["source_commit"] = "0" * 40
        self.assertIn("source_commit", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["object"] = "objects/" + "9" * 64 + ".jsonl.gz"
        self.assertIn("frozen manifest", guard.receipt_error(worker, worker_config, config))
        orbit_config = dict(config, domain_encoding="frobenius_orbit_phase")
        worker_config, worker = receipts(orbit_config, "UNKNOWN_solver_cap")
        self.assertIsNone(guard.receipt_error(worker, worker_config, orbit_config))
        worker_config["strategy"] = "chain-s3"
        self.assertIn("strategy", guard.receipt_error(worker, worker_config, orbit_config))
        ordered_config = dict(config, domain_encoding="frobenius_orbit_phase_ordered_first_pair")
        worker_config, worker = receipts(ordered_config, "UNKNOWN_solver_cap")
        self.assertIsNone(guard.receipt_error(worker, worker_config, ordered_config))
        worker_config["strategy"] = "orbit-chain-s3"
        self.assertIn("strategy", guard.receipt_error(worker, worker_config, ordered_config))

    def test_relation_replay_checks_phase_sign_and_canonical_residue(self):
        self.assertEqual((guard.PRIMARY_FROBENIUS_LAMBDA**2 +
                          guard.PRIMARY_FROBENIUS_LAMBDA + 2) % guard.PRIMARY_ORDER, 0)
        self.assertEqual(pow(guard.PRIMARY_FROBENIUS_LAMBDA, 83, guard.PRIMARY_ORDER), 1)
        config = dict(self.config, max_trials=1)
        worker_config, worker = receipts(config, "PASS_verified_target_only")
        worker["attempt_records"][0]["decomposition_indices"] = [0, 2, 1, 166, 167]
        relation = worker["relation_matrix"][0]
        relation["summands"] = [[0, 0], [0, 1], [0, 0], [1, 0], [1, 0]]
        relation["summand_negated"] = [False, False, True, False, True]
        relation["row"] = ["254512724090651164922414"] + ["0"] * (config["columns"] - 1)
        self.assertIsNone(guard.receipt_error(worker, worker_config, config))
        relation["row"][0] = str(guard.PRIMARY_ORDER + guard.PRIMARY_FROBENIUS_LAMBDA)
        self.assertIn("malformed", guard.receipt_error(worker, worker_config, config))
        relation["row"][0] = "²"
        self.assertIn("malformed", guard.receipt_error(worker, worker_config, config))

    def test_complete_attempt_ledger_requires_nonzero_target_coefficient(self):
        config = dict(self.config, max_trials=1)
        worker_config, worker = receipts(config, "UNKNOWN_solver_cap")
        worker["attempt_records"].clear()
        self.assertIn("incomplete", guard.receipt_error(worker, worker_config, config))
        worker_config, worker = receipts(config, "UNKNOWN_solver_cap")
        worker["attempt_records"][0]["b"] = "0"
        self.assertIn("malformed", guard.receipt_error(worker, worker_config, config))

    def test_independent_group_replay_checks_queries_and_selected_points(self):
        generator = guard.PRIMARY_GENERATOR
        self.assertIsNone(guard.point_mul(generator, guard.PRIMARY_ORDER))
        published = Path(__file__).parent / "pilot-01"
        corpus = json.loads((published / "probe-corpus.json").read_text())
        public_point = tuple(map(int, corpus["targets"][0]["point"]))
        self.assertEqual(guard.point_mul(generator, 1172698111487078683838919), public_point)
        with tempfile.TemporaryDirectory(prefix="n83-group-replay-") as temporary:
            root = Path(temporary)
            fixture = root / "probe-corpus.json"
            target = guard.point_mul(generator, 5)
            fixture.write_text(json.dumps({
                "schema": "n83.public-probe-corpus/v1", "study": guard.STUDY,
                "targets": [{"a": 0, "fixture": 0, "point": list(map(str, target))}],
            }))
            selected = root / "base.jsonl.gz"
            header = {
                "schema": "n83.factor-base/v1", "study": guard.STUDY,
                "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
                "curve_a": 0, "curve_b": "1", "n": 83,
                "modulus_low_terms": [0, 1, 2, 45],
                "closure": "signed_frobenius", "cofactor": "4",
                "orbit_columns": 64, "policy": "public_x_hash", "seed": guard.SEEDS[0],
                "frobenius_lambda": str(guard.PRIMARY_FROBENIUS_LAMBDA),
                "generator": list(map(str, generator)),
                "representatives": [list(map(str, generator))] * 64,
                "point_count": 166 * 64, "subgroup_order": str(guard.PRIMARY_ORDER),
            }
            def write_point(point, include_inverse=False):
                with gzip.open(selected, "wt", encoding="utf-8") as output:
                    output.write(json.dumps(header) + "\n")
                    output.write(json.dumps({
                        "column": 0, "phase": 0, "negative": False,
                        "point": list(map(str, point)),
                    }) + "\n")
                    if include_inverse:
                        output.write(json.dumps({
                            "column": 0, "phase": 0, "negative": True,
                            "point": list(map(str, (generator[0], generator[0] ^ generator[1]))),
                        }) + "\n")
            write_point(generator)
            worker = {"status": "PASS_verified_target_only", "verified_log": "5",
                      "orbit_columns": 64, "policy": "public_x_hash", "seed": guard.SEEDS[0],
                      "report": {"factor_base_size": 166 * 64},
                      "attempt_records": [{
                          "trial": 1, "a": "0", "b": "1",
                          "target": list(map(str, target)),
                          "disposition": "RelationFound",
                          "decomposition_indices": [0] * 5,
                      }]}
            self.assertEqual(guard.replay_group(worker, fixture, selected),
                             {"status": "PASS", "queries": 1, "relations": 1})
            worker["verified_log"] = "6"
            with self.assertRaisesRegex(ValueError, "verified target scalar"):
                guard.replay_group(worker, fixture, selected)
            worker["verified_log"] = "5"
            worker["attempt_records"][0]["target"] = list(map(str, generator))
            with self.assertRaisesRegex(ValueError, "query target"):
                guard.replay_group(worker, fixture, selected)
            worker["attempt_records"][0]["target"] = list(map(str, target))
            write_point(guard.point_mul(generator, 2))
            with self.assertRaisesRegex(ValueError, "signed Frobenius orbit"):
                guard.replay_group(worker, fixture, selected)
            write_point(generator, include_inverse=True)
            worker["attempt_records"][0]["decomposition_indices"] = [0, 0, 0, 0, 1]
            with self.assertRaisesRegex(ValueError, "relation point sum"):
                guard.replay_group(worker, fixture, selected)

    def test_outer_retains_preflight_timeout_resource_and_bad_receipts(self):
        with tempfile.TemporaryDirectory(prefix="n83-chain-search-guard-") as temporary:
            root = Path(temporary)
            panel = root / "panel"
            panel.mkdir()
            (panel / "objects").mkdir()
            bases = []
            for columns, digest, point_set in ((64, "b" * 64, "c" * 64),
                                               (2048, "d" * 64, "e" * 64)):
                name = f"objects/{digest}.jsonl.gz"
                (panel / name).write_bytes(f"synthetic {columns}".encode())
                bases.append({
                    "a": 0, "columns": columns, "policy": "public_x_hash",
                    "seed": guard.SEEDS[0], "object": name,
                    "compressed_blake3": digest, "point_set_blake3": point_set,
                    "s3_uri": f"s3://crypto-autoresearcher/factor-bases/icv1/etc/{guard.STUDY}/a0/{name}",
                })
            (panel / "manifest.json").write_text(json.dumps({
                "schema": "n83.factor-base-panel/v1", "study": guard.STUDY, "bases": bases,
            }))
            (panel / "replay.json").write_text("{}")
            (panel / "upload-receipt.json").write_text("{}")
            fixtures = root / "fixtures"
            fixtures.mkdir()
            (fixtures / "probe-corpus.json").write_text("{}")
            (fixtures / "probe-validation.json").write_text("{}")
            binary = root / "worker"
            binary.write_bytes(b"\x7fELFsynthetic")
            binary.chmod(0o755)
            checkout = Path(guard.__file__).resolve().parents[2]
            expected = {
                "preflight": "PREFLIGHT_ONLY",
                "timeout": "UNKNOWN_wall_cap",
                "resource": "UNKNOWN_resource_or_worker_exit",
                "malformed": "PRODUCER_FAILURE_worker_receipt",
                "source_changed": "PRODUCER_FAILURE_frozen_input_changed",
                "orbit_preflight": "PREFLIGHT_ONLY",
                "orbit_worker_cap": "UNKNOWN_wall_cap",
                "orbit_bad_cap": "PRODUCER_FAILURE_worker_exit",
                "orbit_bad_strategy": "PRODUCER_FAILURE_worker_receipt",
                "ordered_preflight": "PREFLIGHT_ONLY",
                "ordered_bad_strategy": "PRODUCER_FAILURE_worker_receipt",
                "input_changed": "PRODUCER_FAILURE_frozen_input_changed",
                "object_changed": "PRODUCER_FAILURE_frozen_input_changed",
            }
            for mode, status in expected.items():
                with self.subTest(mode=mode):
                    output = root / mode
                    orbit_mode = mode.startswith(("orbit_", "ordered_"))
                    ordered_mode = mode.startswith("ordered_")
                    columns = 2048 if orbit_mode else 64
                    max_variables = 150000 if orbit_mode else 100000
                    max_domain_clauses = 1000000 if orbit_mode else 300000
                    memory_mib = 4096 if orbit_mode else 256

                    def fake_run(command, **_kwargs):
                        if command[:3] == ["docker", "image", "inspect"]:
                            return subprocess.CompletedProcess(command, 0, "sha256:" + "f" * 64 + "\n", "")
                        return subprocess.CompletedProcess(command, 0, "", "")

                    class FakeChild:
                        calls = 0

                        def wait(self, timeout=None):
                            self.calls += 1
                            if mode == "timeout" and self.calls == 1:
                                raise subprocess.TimeoutExpired("docker run", timeout)
                            if mode == "resource":
                                return 137
                            return 124 if mode in ("orbit_worker_cap", "orbit_bad_cap") else 0

                        def kill(self):
                            pass

                    def fake_popen(command, **_kwargs):
                        config = json.loads((output / "config.json").read_text())
                        orbit = orbit_mode
                        self.assertIn("--network", command)
                        self.assertEqual(command[command.index("--network") + 1], "none")
                        self.assertEqual(command[command.index("--memory-swap") + 1], f"{memory_mib}m")
                        selected = ("primary-ordered-orbit-chain-cold" if ordered_mode else
                                    "primary-orbit-chain-cold" if orbit else "primary-chain-cold")
                        self.assertEqual(command[-16:-14], ["/worker", selected])
                        if mode not in ("timeout", "resource"):
                            worker_config, worker = receipts(config)
                            if mode == "malformed":
                                worker["report"]["sat_invalid_models"] = 1
                            if mode in ("orbit_bad_strategy", "ordered_bad_strategy"):
                                worker_config["strategy"] = "chain-s3"
                            run_dir = output / "run"
                            run_dir.mkdir()
                            (run_dir / "config.json").write_text(json.dumps(worker_config))
                            (run_dir / "summary.json").write_text(json.dumps(worker))
                            if mode in ("orbit_worker_cap", "orbit_bad_cap"):
                                cap = {
                                    "schema": "n83.primary-cold-cap/v1",
                                    "status": "UNKNOWN_budget", "orbit_columns": config["columns"],
                                    "summands": config["summands"], "max_trials": config["max_trials"],
                                    "budget_seconds": config["worker_wall_seconds"],
                                    "strategy": "chain-s3" if mode == "orbit_bad_cap" else "orbit-chain-s3",
                                }
                                (run_dir / "cap.json").write_text(json.dumps(cap))
                        if mode == "source_changed":
                            (output / "source" / guard.SOURCES[2]).write_bytes(b"changed")
                        if mode == "input_changed":
                            (fixtures / "probe-corpus.json").write_text('{"changed":true}')
                        if mode == "object_changed":
                            (panel / config["selected_object"]).write_bytes(b"changed object")
                        return FakeChild()

                    with patch.object(guard, "clean_commit", return_value="a" * 40), \
                         patch.object(guard.subprocess, "run", side_effect=fake_run), \
                         patch.object(guard.subprocess, "Popen", side_effect=fake_popen):
                        outer = guard.run(
                            panel, fixtures, columns, "public_x_hash", guard.SEEDS[0], 5,
                            0, max_variables, max_domain_clauses, 1, 1000, 1.0, memory_mib,
                            binary, output, checkout=checkout,
                            orbit_domain=orbit_mode, ordered_first_pair=ordered_mode,
                        )
                    self.assertEqual(outer["status"], status)
                    self.assertEqual(json.loads((output / "outer.json").read_text())["status"], status)


if __name__ == "__main__":
    unittest.main()
