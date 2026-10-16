"""No-launch controls for the paired gate's authorization and receipts."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import blake3

sys.path.insert(0, str(Path(__file__).resolve().parent))
import receipt_audit as audit
import supervisor


class PairedSupervisorControls(unittest.TestCase):
    def setUp(self):
        self.plan = json.loads(supervisor.PLAN.read_text())
        self.arm = self.plan["arms"][0]

    def test_specific_grant_required_before_any_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            auth = root / "auth.json"
            auth.write_text("{}")
            with mock.patch.object(supervisor.PRIOR, "stage") as stage:
                with self.assertRaisesRegex(ValueError, "paired budget schema mismatch"):
                    supervisor.run(auth, root / "out", root, checkout=supervisor.ROOT)
                stage.assert_not_called()
            self.assertFalse((root / "out").exists())

    def test_grant_contract_and_exact_two_cold_arms(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "out"
            plan_hash = supervisor.hashlib.sha256(supervisor.PLAN.read_bytes()).hexdigest()
            record = {
                "schema": "n83.additional-local-budget/v7",
                "study": "koblitz_n83_factor_base_sweep_20261008",
                "purpose": "k600_gray_hash_paired_first_trial",
                "prior_conservative_charged_seconds": "23100",
                "additional_authorized_seconds": "4800",
                "plan_sha256": plan_hash, "output_dir": str(output),
                "user_approval_reference": "explicit later 4800-second grant",
            }
            self.assertIsNone(supervisor.authorization_error(record, output, plan_hash))
            record["plan_sha256"] = "c0f7181b93f7b86bcbbde10d3e9fab10ea82632610d4044ad6a672314b0cc527"
            self.assertEqual(supervisor.authorization_error(record, output, plan_hash),
                             "paired budget plan_sha256 mismatch")
            record["plan_sha256"] = plan_hash
            record["additional_authorized_seconds"] = "2400"
            self.assertEqual(supervisor.authorization_error(record, output, plan_hash),
                             "paired budget additional_authorized_seconds mismatch")
        self.assertEqual(supervisor.ARM_CAP, 2400)
        self.assertEqual(supervisor.TOTAL_CAP, 4800)
        self.assertEqual([arm["role"] for arm in self.plan["arms"]],
                         ["gray_pool", "hash_pool"])
        self.assertEqual([arm["base_object_blake3"] for arm in self.plan["arms"]],
                         [supervisor.preflight.BASE] * 2)

    def test_verify_receipt_rejects_wrong_gray_policy_and_overlap(self):
        record = {
            "schema": "n83.cross-stream-pilot-import/v1", "status": "PASS",
            "curve_slug": self.plan["curve_slug"],
            "small_object_blake3": self.arm["base_object_blake3"],
            "small_policy": "public_x_sequential", "small_seed": 2026100802,
            "small_columns": 600, "small_points": 99600,
            "pool_object_blake3": self.arm["pool_object_blake3"],
            "pool_policy": "public_x_gray_prefix", "pool_seed": 2026100801,
            "pool_columns": 600, "pool_points": 99600,
            "shared_signed_orbits": 503, "small_only_signed_orbits": 97,
            "residual_signed_orbits": 97, "residual_points": 16102,
            "solver_search_executed": False,
        }
        audit.audit_verify(record, self.arm)
        record["pool_policy"] = "public_x_gray"
        with self.assertRaisesRegex(ValueError, "verify pool_policy changed"):
            audit.audit_verify(record, self.arm)
        record["pool_policy"] = "public_x_gray_prefix"
        record["shared_signed_orbits"] = 502
        with self.assertRaisesRegex(ValueError, "verify shared_signed_orbits changed"):
            audit.audit_verify(record, self.arm)

    def test_worker_commands_preserve_paired_domains_and_isolation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = self.plan["arms"][0]
            gray, _, gray_name = supervisor.stage_command(
                "search", base, supervisor.ROOT, root, root / "gray_build",
                root / "gray_run", root / "gray_worker", root / "registry",
                supervisor.PRIOR.IMAGE, "fixedtoken")
            other = self.plan["arms"][1]
            hashed, _, hash_name = supervisor.stage_command(
                "search", other, supervisor.ROOT, root, root / "hash_build",
                root / "hash_run", root / "hash_worker", root / "registry",
                supervisor.PRIOR.IMAGE, "fixedtoken")
        self.assertNotEqual(gray_name, hash_name)
        self.assertIn("--network", gray)
        self.assertIn("none", gray)
        self.assertIn("--cpus", gray)
        self.assertIn("1", gray)
        self.assertIn("/gray_worker", " ".join(gray))
        self.assertIn("/hash_worker", " ".join(hashed))
        self.assertEqual(gray[-len(base["worker_search_argv"]):], base["worker_search_argv"])
        self.assertEqual(hashed[-len(other["worker_search_argv"]):], other["worker_search_argv"])

    def test_search_receipt_checks_source_and_censored_status(self):
        manifest = supervisor.preflight.PINNED_SHA256["manifest"][0].read_bytes()
        source = supervisor.preflight.PINNED_SHA256["worker_source"][0].read_bytes()
        config = {
            "schema": "n83.cross-stream-pilot-partial-config/v1",
            "study": self.plan["study"],
            "base_object_blake3": self.arm["base_object_blake3"],
            "pool_object_blake3": self.arm["pool_object_blake3"],
            "pilot_manifest_blake3": blake3.blake3(manifest).hexdigest(),
            "source_blake3": blake3.blake3(source).hexdigest(),
            "public_corpus_path": "/corpus.json",
            "public_corpus_canonical_json_blake3": audit.PRIOR.CORPUS_BLAKE3,
            "budget_seconds": 1100, "conflict_budget": 5000,
            "max_models": 32, "base_summands": 3, "residual_summands": 2,
            "order_first_pair": True, "max_variables": 150000,
            "max_domain_clauses": 12000000,
            "memory_cgroup_limit_bytes": 4294967296,
            "memory_cgroup_swap_limit_bytes": 0,
        }
        result = {
            "schema": "n83.cross-stream-pilot-partial-result/v1",
            "study": self.plan["study"], "status": "UNKNOWN_search_cap",
            "base_object_blake3": self.arm["base_object_blake3"],
            "pool_object_blake3": self.arm["pool_object_blake3"],
            "base_policy": "public_x_sequential", "base_seed": 2026100802,
            "pool_policy": "public_x_gray_prefix", "pool_seed": 2026100801,
            "public_corpus_canonical_json_blake3": audit.PRIOR.CORPUS_BLAKE3,
            "base_signed_orbits": 600, "pool_signed_orbits": 600,
            "shared_signed_orbits": 503, "residual_signed_orbits": 97,
            "base_import_elapsed_ms_L0": 20.0, "worker_wall_ms_L0": 40.0,
            "memory_cgroup_peak_bytes": 50000000,
            "trial": 0, "a": self.plan["first_trial_a"],
            "b": self.plan["first_trial_b"],
            "search": {"domain_clauses": self.arm["source_capacity"]["domain_clauses"],
                       "peak_variables": 1000, "peak_xors": 100,
                       "peak_clauses": 2000, "conflicts": 5000,
                       "models": 0, "invalid_models": 0,
                       "group_rejected_models": 0, "solver_calls": 1,
                       "patterns_attempted": 1, "exhausted": True,
                       "refuted": False},
            "candidate": None, "graph_rows": None, "matrix_rank": None,
            "total_index_calculus_runtime_ms": None,
            "selected_best_total_runtime": None,
        }
        self.assertEqual(audit.audit_search(config, result, self.arm, [], set()),
                         "UNKNOWN_search_cap")
        config["source_blake3"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "search config source_blake3 changed"):
            audit.audit_search(config, result, self.arm, [], set())
        config["source_blake3"] = blake3.blake3(source).hexdigest()
        result["status"] = "PASS_one_trial_refuted"
        with self.assertRaisesRegex(ValueError, "refutation counter changed"):
            audit.audit_search(config, result, self.arm, [], set())


if __name__ == "__main__":
    unittest.main()
