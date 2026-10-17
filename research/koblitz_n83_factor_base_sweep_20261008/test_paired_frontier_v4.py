"""Exact ordinal, source, and fail-closed controls for the N83 v4 grid."""

import json
import unittest
from unittest.mock import patch

import paired_frontier_v4 as gate


class PairedFrontierV4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent = json.loads(gate.PARENT.read_text())
        cls.design = json.loads(gate.OUTPUT.read_text())
        cls.base = json.loads(gate.v3.PARENT.read_text())
        cls.keys = gate.pilot_keys(json.loads(gate.PILOT.read_text()))

    def test_saved_population_and_source_classes(self):
        self.assertEqual(gate.v3.canonical_bytes(gate.build()), gate.OUTPUT.read_bytes())
        self.assertEqual(self.design["stream_count"], 54)
        self.assertEqual(self.design["cross_stream_pair_count"], 169884)
        self.assertEqual(self.design["cross_stream_tuple_count"], 26420359680)
        self.assertEqual(self.design["total_tuple_count"], 26939278080)
        self.assertEqual(sum(self.design["source_class_pair_counts"].values()), 169884)
        self.assertEqual(sum(self.design["disposition_counts"].values()), 26420359680)
        self.assertEqual(self.design["source_class_pair_counts"][
            "pilot_exact_import_gate_available_unexecuted"], 648)
        self.assertEqual(self.design["disposition_counts"][
            "pilot_exact_import_gate_available_unexecuted"], 2592)

    def test_v3_prefix_is_byte_equivalent(self):
        for ordinal in (0, 57023999, 57024000, 518918399):
            self.assertEqual(gate.case(self.design, self.parent, ordinal),
                             gate.v3.case(self.parent, self.base, ordinal))

    def test_every_pair_resolves_to_distinct_frozen_streams(self):
        counts = {}
        for index in range(gate.PAIR_COUNT):
            spec = gate.pair(self.design, index)
            self.assertEqual(spec["curve_a"], self.base["base_specs"][
                spec["base"]["parent_base_index"]]["a"])
            self.assertEqual(spec["curve_a"], self.base["base_specs"][
                spec["pool"]["parent_base_index"]]["a"])
            for name in ("base", "pool"):
                member = spec[name]
                frozen = self.base["base_specs"][member["parent_base_index"]]
                self.assertEqual((member["policy"], member["seed"],
                                  member["closure"], member["orbits"]),
                                 (frozen["family"], frozen["seed"],
                                  frozen["closure"], frozen["parameter"]))
            self.assertNotEqual(
                (spec["base"]["policy"], spec["base"]["seed"],
                 spec["base"]["closure"]),
                (spec["pool"]["policy"], spec["pool"]["seed"],
                 spec["pool"]["closure"]),
            )
            category = gate.source_class(spec, self.keys)
            counts[category] = counts.get(category, 0) + 1
        self.assertEqual(counts, self.design["source_class_pair_counts"])

    def test_boundary_and_pilot_case_leave_overlap_unmeasured(self):
        first = gate.case(self.design, self.parent, 518918400)
        last = gate.case(self.design, self.parent, 26939278079)
        self.assertEqual((first["pair_index"], first["residual_index"],
                          first["solver_index"]), (0, 0, 0))
        self.assertEqual((last["pair_index"], last["residual_index"],
                          last["solver_index"]), (169883, 2, 51839))
        with self.assertRaisesRegex(ValueError, "outside v4 grid"):
            gate.case(self.design, self.parent, 26939278080)
        first_pilot = next(i for i in range(gate.PAIR_COUNT)
                           if gate.source_class(gate.pair(self.design, i), self.keys)
                           == "pilot_exact_import_gate_available_unexecuted")
        solver = {"summands": 5, "split_bits": 0,
                  "symmetry": gate.v3.EXACT_PAIR_SYMMETRY,
                  "backend": "sat_cdcl", "enumeration": "binary",
                  "max_large_primes": 2, "threads": 1,
                  "linear_algebra": "modular_biguint"}
        solver_index = 0
        for key in gate.v3.v2lib.AXIS_KEYS:
            axis = self.parent["solver_axes"][key]
            solver_index = solver_index * len(axis) + axis.index(solver[key])
        ordinal = (self.design["v3_tuple_count_preserved"]
                   + first_pilot * 3 * 51840 + 2 * 51840 + solver_index)
        result = gate.case(self.design, self.parent, ordinal)
        self.assertEqual(result["disposition"],
                         "pilot_exact_import_gate_available_unexecuted")
        self.assertFalse(result["exact_overlap_measured"])
        self.assertIsNone(result["conditional_unordered_point_multisets"])
        self.assertIsNone(result["uniform_nonidentity_hit_ceiling"])
        self.assertIsNone(result["selected_best_total_runtime"])
        solver["backend"] = "wdsat"
        self.assertEqual(gate.disposition(result["pair"], solver, 2, self.keys),
                         "paired_backend_adapter_unavailable:wdsat")

    def test_frozen_parent_and_pilot_are_required(self):
        with patch.object(gate, "PARENT_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "v3 parent bytes changed"):
                gate.build()
        with patch.object(gate, "PILOT_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "pilot source census bytes changed"):
                gate.build()


if __name__ == "__main__":
    unittest.main()
