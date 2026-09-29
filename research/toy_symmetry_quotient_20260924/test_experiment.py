import itertools
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest

from experiment import (O, ToyCurve, candidates, catalogue, digest, fixtures,
                        frobenius_audit, lift_roots, permutation_audit,
                        require_closed, root_polynomial, verify_catalogue,
                        write_new_results)


class SymmetryPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.curve = ToyCurve()
        cls.points = cls.curve.points()
        cls.ids = {p: i for i, p in enumerate(cls.points)}
        cls.cases = {c["id"]: c for c in fixtures(cls.curve, cls.points)}

    def test_field_inverses_and_distributivity(self):
        c = ToyCurve()
        for a in range(1, 128):
            self.assertEqual(c.mul(a, c.inv(a)), 1)
        rng = random.Random(19)
        for _ in range(256):
            a, b, d = [rng.randrange(128) for _ in range(3)]
            self.assertEqual(c.mul(a, b ^ d), c.mul(a, b) ^ c.mul(a, d))
        with self.assertRaises(ZeroDivisionError):
            c.inv(0)

    def test_group_order_against_frobenius_trace_recurrence(self):
        prev, current = 2, -1
        for _ in range(2, 8):
            prev, current = current, -current - 2 * prev
        self.assertEqual(len(self.points), 128 + 1 - current)

    def test_every_group_addition_against_cyclic_table(self):
        c = ToyCurve()
        cycle = None
        for generator in self.points[1:]:
            sequence, seen, p = [], set(), O
            while p not in seen:
                seen.add(p)
                sequence.append(p)
                p = c.add(p, generator)
            if len(sequence) == len(self.points):
                cycle = sequence
                break
        self.assertIsNotNone(cycle)
        self.assertEqual(set(cycle), set(self.points))
        for i, p in enumerate(cycle):
            for j, q in enumerate(cycle):
                self.assertEqual(c.add(p, q), cycle[(i + j) % len(cycle)])

    def test_frobenius_is_group_automorphism(self):
        c = ToyCurve()
        for p in self.points:
            actual = p
            for _ in range(7):
                actual = c.frobenius(actual)
            self.assertEqual(actual, p)
            self.assertTrue(c.contains(c.frobenius(p)))
        for p, q in itertools.product(self.points, repeat=2):
            self.assertEqual(c.frobenius(c.add(p, q)), c.add(c.frobenius(p), c.frobenius(q)))

    def test_subspace_dimension_does_not_predict_point_yield(self):
        c = ToyCurve()
        for middle, expected_points in ((2, 1), (4, 15)):
            case = self.cases[f"linear-d3-x{middle}"]
            v = set(case["subspace"])
            self.assertEqual(len(v), 8)
            self.assertEqual({a ^ b for a in v for b in v}, v)
            self.assertEqual({c.square(x) for x in v}, v)
            self.assertEqual(len(case["support"]), expected_points)
            require_closed(c, case["support"])

    def test_noninvariant_support_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not Frobenius invariant"):
            require_closed(ToyCurve(), self.cases["ordinary-s101"]["support"])

    def test_exact_coverage_on_frozen_and_holdout_supports(self):
        for name in ("orbit-s101", "ordinary-s503"):
            support = self.cases[name]["support"]
            reference = catalogue(ToyCurve(), support, "ordered")
            for variant in ("permutation", "invariant_lift"):
                c = ToyCurve()
                actual = catalogue(c, support, variant)
                self.assertEqual(actual, reference)
                verify_catalogue(c, actual, support)
            audit = permutation_audit(support, reference)
            self.assertEqual(audit["unordered"], 680)
            self.assertEqual(audit["reconstructed_ordered"], 3375)

    def test_invariants_preserve_multiplicity_and_permutations(self):
        c = ToyCurve()
        for roots in ((0, 0, 0), (1, 1, 3), (2, 5, 7)):
            poly = root_polynomial(c, roots)
            for order in itertools.permutations(roots):
                self.assertEqual(root_polynomial(c, order), poly)
            self.assertEqual(lift_roots(c, poly, range(128)), tuple(sorted(roots)))

    def test_invariant_lifting_rejects_missing_roots(self):
        c = ToyCurve()
        poly = root_polynomial(c, (1, 2, 3))
        with self.assertRaisesRegex(ValueError, "does not split"):
            lift_roots(c, poly, (1, 2))
        with self.assertRaisesRegex(ValueError, "monic cubic"):
            lift_roots(c, (1, 1, 0, 0), (1,))

    def test_point_lifting_does_not_forget_y(self):
        c = ToyCurve()
        p, q = (1, 0), (1, 1)
        support = (p, q)
        lifts = list(candidates(c, support, "invariant_lift"))
        self.assertEqual(set(lifts), set(itertools.product(support, repeat=3)))
        self.assertNotEqual(c.total((p, p, p)), c.total((q, p, p)))

    def test_whole_target_transport_and_burnside(self):
        support = self.cases["orbit-s101"]["support"]
        truth = catalogue(ToyCurve(), support, "ordered")
        audit = frobenius_audit(ToyCurve(), truth, support, self.points, self.ids)
        self.assertEqual(audit["status"], "VERIFIED")
        self.assertEqual(audit["target_orbits"], 20)
        bad = audit["independent_rotation_counterexample"]
        self.assertNotEqual(bad["target"], bad["changed_target"])
        self.assertEqual(sum(audit["fixed_solution_pairs_by_power"]), 7 * audit["solution_orbits"])

    def test_singleton_support_retains_unsatisfiable_targets(self):
        support = ((0, 1),)
        truth = catalogue(ToyCurve(), support, "permutation")
        self.assertEqual(len(truth), 1)
        self.assertEqual(len(self.points) - len(truth), 115)
        audit = frobenius_audit(ToyCurve(), truth, support, self.points, self.ids)
        self.assertEqual(audit["negative_control_status"], "TRIVIAL_ACTION")

    def test_bad_decomposition_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, "invalid point decomposition"):
            verify_catalogue(ToyCurve(), {O: {((0, 1),) * 3}}, ((0, 1),))

    def test_support_bounds_and_domain(self):
        for support in ((), (O,), ((0, 0),), ((0, 1), (0, 1)), self.points[1:26]):
            with self.assertRaises(ValueError):
                catalogue(ToyCurve(), support, "ordered")

    def test_fixtures_reproduce(self):
        self.assertEqual(digest(list(self.cases.values())), digest(fixtures(ToyCurve(), self.points)))
        frozen = {tuple(c["support"]) for c in self.cases.values() if c["split"] == "frozen"}
        holdout = {tuple(c["support"]) for c in self.cases.values() if c["split"] == "holdout"}
        self.assertFalse(frozen & holdout)

    def test_existing_results_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileExistsError):
                write_new_results(Path(directory), 1)

    def test_committed_evidence_integrity(self):
        root = Path(__file__).parent / "results" / "pilot-001"
        report = json.loads((root / "results.json").read_text())
        instances = json.loads((root / "instances.json").read_text())
        certs = json.loads((root / "certificates.json").read_text())
        self.assertEqual(digest(instances), report["instance_sha256"])
        self.assertEqual(digest(certs), report["certificate_sha256"])
        self.assertEqual(hashlib.sha256((root / "source.py").read_bytes()).hexdigest(), report["source_sha256"])
        self.assertEqual(len(report["measurements"]), 90)
        for row in report["measurements"]:
            self.assertEqual(row["certificate_sha256"], certs[row["case"]]["sha256"])
            self.assertEqual(row["status"], "VERIFIED_COMPLETE_ENUMERATION")
        for cert in certs.values():
            self.assertEqual(digest(cert["catalogue"]), cert["sha256"])
        for row in report["summary"]:
            self.assertTrue(row["correct"])
            self.assertIsNone(row["max_processed_degree"])
            self.assertIsNone(row["full_dlp_speedup"])
            self.assertIsNone(row["ratio_to_rho"])


if __name__ == "__main__":
    unittest.main()
