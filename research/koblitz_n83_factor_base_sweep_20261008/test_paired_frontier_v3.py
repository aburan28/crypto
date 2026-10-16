"""Exact addressability and provenance controls for the paired N83 grid."""

import json
import unittest
from unittest.mock import patch

import paired_frontier_v3 as gate


class PairedFrontierV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent = json.loads(gate.PARENT.read_text())
        cls.design = json.loads(gate.OUTPUT.read_text())

    def ordinal(self, pair_index: int, residual: int, solver: dict) -> int:
        solver_index = 0
        for key in gate.v2lib.AXIS_KEYS:
            values = self.design["solver_axes"][key]
            solver_index = solver_index * len(values) + values.index(solver[key])
        return (self.design["v2_tuple_count_preserved"]
                + pair_index * 3 * self.design["solver_tuple_count"]
                + residual * self.design["solver_tuple_count"] + solver_index)

    def test_saved_design_and_exact_population(self):
        self.assertEqual(gate.canonical_bytes(gate.build()), gate.OUTPUT.read_bytes())
        self.assertEqual(self.design["pair_spec_count"], 2970)
        self.assertEqual(self.design["solver_tuple_count"], 51840)
        self.assertEqual(self.design["paired_tuple_count"], 461894400)
        self.assertEqual(self.design["total_tuple_count"], 518918400)
        self.assertEqual(sum(self.design["disposition_counts"].values()), 461894400)
        self.assertEqual(self.design["disposition_counts"][
            "bounded_mixed_orbit_partial_worker_available_full_pipeline_unmeasured"], 8)

    def test_v2_ordinals_are_exact_prefix(self):
        for ordinal in (0, 43199, 43200, 45359999, 57023999):
            self.assertEqual(gate.case(self.design, self.parent, ordinal),
                             gate.v2lib.case(self.parent, ordinal))

    def test_pair_and_solver_boundary_bijection(self):
        first = gate.case(self.design, self.parent, 57024000)
        self.assertEqual((first["pair_index"], first["residual_index"],
                          first["solver_index"]), (0, 0, 0))
        last = gate.case(self.design, self.parent, 518918399)
        self.assertEqual((last["pair_index"], last["residual_index"],
                          last["solver_index"]), (2969, 2, 51839))
        self.assertEqual(first["disposition"], "no_residual_pair_axis_redundant")
        with self.assertRaisesRegex(ValueError, "outside v3 grid"):
            gate.case(self.design, self.parent, 518918400)
        with self.assertRaisesRegex(ValueError, "outside v3 grid"):
            gate.case(self.design, self.parent, -1)

    def test_every_pair_resolves_to_the_frozen_parent_bases(self):
        seen = set()
        for pair in self.design["pair_specs"]:
            small = self.parent["base_specs"][pair["small_base_index"]]
            envelope = self.parent["base_specs"][pair["envelope_base_index"]]
            for key, base_key in (("curve_a", "a"), ("policy", "family"),
                                  ("seed", "seed"), ("closure", "closure")):
                self.assertEqual(pair[key], small[base_key])
                self.assertEqual(pair[key], envelope[base_key])
            self.assertEqual(pair["small_orbits"], small["parameter"])
            self.assertEqual(pair["envelope_orbits"], envelope["parameter"])
            self.assertLess(pair["small_orbits"], pair["envelope_orbits"])
            identity = (pair["small_base_index"], pair["envelope_base_index"])
            self.assertNotIn(identity, seen)
            seen.add(identity)
        self.assertEqual(len(seen), 2970)

    def test_exact_worker_pairs_and_conditional_support(self):
        available = [(i, pair) for i, pair in enumerate(self.design["pair_specs"])
                     if gate.pair_class(pair) == "source_available"]
        self.assertEqual([(p["small_orbits"], p["envelope_orbits"])
                          for _, p in available], [(600, 2048), (1182, 2048)])
        solver = {
            "summands": 5, "split_bits": 0,
            "symmetry": gate.EXACT_PAIR_SYMMETRY,
            "backend": "sat_cdcl", "enumeration": "binary",
            "max_large_primes": 2, "threads": 1,
            "linear_algebra": "modular_biguint",
        }
        expected = {600: "4757348390396427029827200",
                    1182: "13009394065076279872900344"}
        for index, pair in available:
            row = gate.case(self.design, self.parent,
                            self.ordinal(index, 2, solver))
            self.assertEqual(row["disposition"],
                             "bounded_mixed_orbit_partial_worker_available_full_pipeline_unmeasured")
            self.assertEqual(row["conditional_unordered_point_multisets"],
                             expected[pair["small_orbits"]])
            self.assertFalse(row["executed_complete_pipeline"])
            self.assertIsNone(row["selected_best_total_runtime"])
        for summands in (5, 6):
            for residual in (1, 2):
                solver["summands"] = summands
                solver["max_large_primes"] = residual
                for index, _ in available:
                    row = gate.case(self.design, self.parent,
                                    self.ordinal(index, residual, solver))
                    self.assertIn("partial_worker_available", row["disposition"])
        solver["backend"] = "wdsat"
        row = gate.case(self.design, self.parent,
                        self.ordinal(available[0][0], 2, solver))
        self.assertEqual(row["disposition"], "paired_backend_adapter_unavailable:wdsat")

    def test_unverified_closure_does_not_inherit_signed_point_count(self):
        index = next(i for i, pair in enumerate(self.design["pair_specs"])
                     if pair["closure"] == "none")
        row = gate.case(self.design, self.parent,
                        57024000 + index * 3 * 51840 + 2 * 51840)
        self.assertIsNone(row["conditional_unordered_point_multisets"])
        self.assertIsNone(row["uniform_nonidentity_hit_ceiling"])

    def test_frozen_parent_and_import_evidence_are_required(self):
        with patch.object(gate, "PARENT_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "parent bytes changed"):
                gate.build()
        with patch.object(gate, "IMPORT_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "import receipt changed"):
                gate.build()
        with patch.object(gate, "SUPPORT_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "support receipt changed"):
                gate.build()


if __name__ == "__main__":
    unittest.main()
