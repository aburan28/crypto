"""Source-only fail-closed controls for the first cross-stream N83 gate."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import supervisor as gate


class CrossStreamFirstGateTests(unittest.TestCase):
    def authorization(self, output: Path) -> dict:
        return {
            "schema": "n83.additional-local-budget/v6",
            "study": gate.STUDY,
            "purpose": "one_k600_hash_vs_sequential_cross_stream_partial",
            "prior_total_charged_seconds": "17400",
            "prior_k600_audit_sha256": gate.reuse.AUDIT_SHA256,
            "base_object_blake3": gate.BASE_HASH,
            "pool_object_blake3": gate.POOL_HASH,
            "additional_authorized_seconds": "2400",
            "output_dir": str(output.resolve()),
            "user_approval_reference": "explicit future approval fixture",
        }

    def test_fresh_pair_specific_grant_is_required(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            good = self.authorization(output)
            self.assertIsNone(gate.authorization_error(good, output))
            for key, value in (
                ("schema", "n83.additional-local-budget/v5"),
                ("prior_total_charged_seconds", "15000"),
                ("prior_k600_audit_sha256", "0" * 64),
                ("base_object_blake3", gate.POOL_HASH),
                ("pool_object_blake3", gate.BASE_HASH),
                ("additional_authorized_seconds", "2399"),
                ("output_dir", str(Path(folder) / "other")),
                ("user_approval_reference", ""),
            ):
                self.assertIsNotNone(gate.authorization_error({**good, key: value}, output))

    def test_stale_grant_cannot_create_output_or_launch(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output, base, pool = (root / item for item in ("run", "base", "pool"))
            base.write_bytes(b"base")
            pool.write_bytes(b"pool")
            stale = self.authorization(output)
            stale["prior_total_charged_seconds"] = "15000"
            auth = root / "auth.json"
            auth.write_text(json.dumps(stale))
            with mock.patch.object(gate.prior, "stage", side_effect=AssertionError("launched")):
                with self.assertRaisesRegex(ValueError, "authorization"):
                    gate.run(auth, output, base, pool)
            self.assertFalse(output.exists())

    def test_pinned_pair_and_prior_ledger(self):
        base, pool = gate.selected_rows(gate.ROOT)
        self.assertEqual(base["compressed_blake3"], gate.BASE_HASH)
        self.assertEqual(pool["compressed_blake3"], gate.POOL_HASH)
        self.assertIsNone(gate.reuse.prior_evidence_error(gate.ROOT))
        self.assertEqual(gate.PRIOR_TOTAL + gate.TOTAL_CHARGE, 19800)
        self.assertEqual(gate.object_error(Path(__file__), Path(__file__), gate.ROOT),
                         "retained K600 base compressed bytes differ")

    def test_orbit_set_difference_rejects_partial_and_empty(self):
        base = [(1, 10), (1, 11), (2, 20), (2, 21)]
        pool = [(2, 20), (2, 21), (3, 30), (3, 31)]
        self.assertEqual(gate.orbit_difference(base, pool, 2),
                         (1, [(3, 30), (3, 31)]))
        with self.assertRaisesRegex(ValueError, "partial signed-orbit"):
            gate.orbit_difference(base, [(2, 20), (3, 31)], 2)
        with self.assertRaisesRegex(ValueError, "residual is empty"):
            gate.orbit_difference(base, base, 2)
        with self.assertRaisesRegex(ValueError, "x-domain"):
            gate.orbit_difference(base, [(1, 30), (1, 31)], 2)

    def test_verify_and_one_trial_receipts_fail_closed(self):
        verify = {
            "schema": "n83.cross-stream-pilot-import/v1", "status": "PASS",
            "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
            "small_object_blake3": gate.BASE_HASH,
            "pool_object_blake3": gate.POOL_HASH,
            "small_policy": "public_x_hash", "small_seed": 2026100801,
            "small_columns": 600, "small_points": 99600,
            "pool_policy": "public_x_sequential", "pool_seed": 2026100801,
            "pool_columns": 600, "pool_points": 99600,
            "shared_signed_orbits": 598, "small_only_signed_orbits": 2,
            "residual_signed_orbits": 2, "residual_points": 332,
            "solver_search_executed": False,
        }
        self.assertIsNone(gate.verify_error(verify, 598, 2))
        self.assertIsNotNone(gate.verify_error({**verify, "residual_signed_orbits": 3}, 598, 2))
        required = gate.capacity.requirements(600, 602, 5, 2)
        config = {
            "schema": "n83.cross-stream-pilot-partial-config/v1",
            "base_object_blake3": gate.BASE_HASH, "pool_object_blake3": gate.POOL_HASH,
            "public_corpus_canonical_json_blake3": gate.prior.CORPUS_BLAKE3,
            "conflict_budget": 5000, "max_models": 32,
            "base_summands": 3, "residual_summands": 2,
            "budget_seconds": gate.WORKER_WALL, "order_first_pair": True,
            "max_variables": 150000, "max_domain_clauses": 12000000,
            "memory_cgroup_limit_bytes": 4294967296,
            "memory_cgroup_swap_limit_bytes": 0,
        }
        result = {
            "schema": "n83.cross-stream-pilot-partial-result/v1",
            "status": "UNKNOWN_search_cap",
            "base_object_blake3": gate.BASE_HASH, "pool_object_blake3": gate.POOL_HASH,
            "base_policy": "public_x_hash", "base_seed": 2026100801,
            "pool_policy": "public_x_sequential", "pool_seed": 2026100801,
            "public_corpus_canonical_json_blake3": gate.prior.CORPUS_BLAKE3,
            "base_signed_orbits": 600, "pool_signed_orbits": 600,
            "shared_signed_orbits": 598, "residual_signed_orbits": 2,
            "trial": 0, "a": gate.prior.FIRST_TRIAL_A,
            "b": gate.prior.FIRST_TRIAL_B,
            "graph_rows": None, "matrix_rank": None,
            "total_index_calculus_runtime_ms": None,
            "selected_best_total_runtime": None,
            "base_import_elapsed_ms_L0": 10.0, "worker_wall_ms_L0": 100.0,
            "memory_cgroup_peak_bytes": 1000,
            "candidate": None,
            "search": {"conflicts": 5000, "models": 0, "invalid_models": 0,
                       "group_rejected_models": 0,
                       "solver_calls": 1, "patterns_attempted": 1,
                       "exhausted": True, "refuted": False,
                       "domain_clauses": required["domain_clauses"],
                       "peak_variables": required["total_variables"],
                       "peak_xors": required["xor_rows"],
                       "peak_clauses": required["domain_clauses"]},
        }
        residual = [(i + 10, i + 100) for i in range(332)]
        self.assertIsNone(gate.search_error(config, result, 598, [], residual, (1, 1)))
        changed = json.loads(json.dumps(result))
        changed["search"]["conflicts"] = 5001
        self.assertIn("conflicts", gate.search_error(config, changed, 598, [], residual, (1, 1)))
        changed = json.loads(json.dumps(result))
        changed["search"]["domain_clauses"] += 1
        self.assertIn("capacity", gate.search_error(config, changed, 598, [], residual, (1, 1)))
        changed = {**result, "status": "PASS_one_bounded_partial"}
        self.assertIn("candidate/status", gate.search_error(config, changed, 598, [], residual, (1, 1)))


if __name__ == "__main__":
    unittest.main()
