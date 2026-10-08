#!/usr/bin/env python3
"""Archived-base regression against the prior independent brute enumerator."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
M3 = ROOT / "research/notes/ecc2k130/m3_four_policy_20260930"
sys.path.insert(0, str(M3))
import verify as prior_verify  # noqa: E402
from score import choose, score_candidate  # noqa: E402


class ScoreRegression(unittest.TestCase):
    def test_locked_entrypoints_import_without_opening_fixtures(self) -> None:
        for filename in ("produce.py", "verify.py"):
            with self.subTest(filename=filename):
                spec = importlib.util.spec_from_file_location(
                    f"selector_{filename[:-3]}", HERE / filename)
                assert spec is not None and spec.loader is not None
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                self.assertEqual(module.checked_config()["base_size"], 8)

    @classmethod
    def setUpClass(cls) -> None:
        prior_verify.checked_config()
        lock = json.loads((M3 / "FROZEN_REPLAY.json").read_text())
        raw = M3 / "evidence_run_36722040881/result.json"
        assert hashlib.sha256(raw.read_bytes()).hexdigest() == lock["result_sha256"]
        cls.archived = json.loads(raw.read_text())
        prior_verify.restore_bare_curve()
        field = prior_verify.pilot.FastGF2m(21, prior_verify.pilot.IRR)
        cls.reference = prior_verify.pilot.Koblitz(field, 0, 1)
        generator = tuple(cls.archived["challenge"]["G"])
        _, by_point, _ = prior_verify.independent_orbits(cls.reference, generator)
        cls.universe = set(by_point)
        assert len(cls.universe) == 420

    def test_archived_bases_match_independent_brute_rows(self) -> None:
        for seed in sorted(self.archived["bases"]):
            for policy in ("original", "pullback"):
                with self.subTest(seed=seed, policy=policy):
                    base = [tuple(point) for point in self.archived[
                        "bases"][seed][policy]]
                    meter = prior_verify.pilot.Meter()
                    field = prior_verify.pilot.CountedField(meter)
                    curve = prior_verify.pilot.CountedCurve(meter, field, 0, 1)
                    before = meter.snapshot()
                    actual = score_candidate(curve, base, self.universe, meter,
                        prior_verify.pilot.RankTracker)
                    outer = meter.delta(before, meter.snapshot())
                    self.assertEqual(actual["score_cost"], {k: v for k, v in
                        outer.items() if k != "cpu_ns"})
                    brute, _ = prior_verify.brute_triples(self.reference, base)
                    brute.pop(None, None)
                    rows = []
                    for target in sorted(brute):
                        k, i, j = brute[target][0]
                        row = [0] * 8
                        for index in (i, j, k):
                            row[index] += 1
                        rows.append((row, 0))
                    rank, _ = prior_verify.independent_linear_system(rows, 8)
                    self.assertEqual(actual["distinct_support"], len(brute))
                    self.assertEqual(actual["first_witness_base_row_rank"], rank)
                    self.assertEqual(actual["distinct_first_witness_rows"],
                                     len({tuple(row) for row, _ in rows}))
                    self.assertEqual(actual["score_cost"]["group_add"], 324)
                    self.assertEqual(len(actual["first_witness_sha256"]), 64)

    def test_selection_rule_is_deterministic(self) -> None:
        first = {"distinct_support": 95,
                 "first_witness_base_row_rank": 8,
                 "distinct_first_witness_rows": 40}
        second = {"distinct_support": 96,
                  "first_witness_base_row_rank": 8,
                  "distinct_first_witness_rows": 39}
        self.assertEqual(choose([first, second]), 1)
        second["first_witness_base_row_rank"] = 7
        self.assertEqual(choose([first, second]), 0)
        first["first_witness_base_row_rank"] = 7
        self.assertIsNone(choose([first, second]))

    def test_replay_cell_on_archived_panel(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "selector_archived_replay", HERE / "verify.py")
        assert spec is not None and spec.loader is not None
        verifier = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(verifier)
        seed = sorted(self.archived["bases"])[0]
        label, policy = "A", "original"
        G = tuple(self.archived["challenge"]["G"])
        Q = tuple(self.archived["challenge"]["Q"])
        base = [tuple(P) for P in self.archived["bases"][seed][policy]]
        targets = [(u, v, self.reference.add(
            self.reference.mul(G, u), self.reference.mul(Q, v)))
            for u, v in self.archived["target_streams"][label]["coefficients"]]
        path = M3 / "evidence_run_36722040881/cells" / seed / label / f"{policy}.json"
        cell = verifier.replay_cell(path,
            self.archived["cell_sha256"][seed][label][policy],
            self.archived["variants"][seed][label][policy],
            self.reference, prior_verify.ReferenceCurve(1), G, Q, targets,
            base, self.archived["challenge"]["secret_audit_only"])
        self.assertEqual(len(cell["cases"]), 512)
        self.assertTrue(cell["verified"])


if __name__ == "__main__":
    unittest.main()
