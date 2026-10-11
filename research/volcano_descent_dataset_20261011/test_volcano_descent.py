#!/usr/bin/env python3
"""Fail-closed tests for the volcano descent library, dataset builder and rho solver.

Run:  python3 -I -m unittest -v test_volcano_descent      (about one minute)
"""
import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from volcano_descent import (  # noqa: E402
    class_number, curve_for_j_fast, decompose, descend_to_floor, expected_level_counts,
    explore_from_crater, find_instance, group_structure, level_counts, pari, predicted_group, v_ell,
)
from ecdlp_rho import Curve, pick_base_point, solve_ph  # noqa: E402

CASES = [  # (ell, D, c, k, m, bits)
    (2, -7, 1, 6, 0, 40), (2, -7, 1, 6, 6, 40), (2, -4, 1, 5, 5, 40), (2, -3, 1, 5, 3, 40),
    (2, -15, 1, 5, 2, 44), (2, -7, 3, 4, 2, 44), (3, -8, 1, 4, 2, 40), (3, -3, 1, 4, 2, 40),
    (5, -4, 1, 3, 1, 40), (7, -11, 1, 3, 3, 40), (3, -8, 2, 3, 1, 44),
]


class VolcanoTests(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(2026)

    def instance(self, case):
        ell, D, c, k, m, bits = case
        return find_instance(ell, D, c, k, m, bits, self.rng)

    def test_instance_algebra(self):
        for case in CASES:
            inst = self.instance(case)
            D, f = decompose(inst.p, inst.t)
            self.assertEqual((D, f), (case[1], case[2] * case[0] ** case[3]))
            self.assertTrue(pari.isprime(inst.p))
            self.assertEqual(v_ell(f, case[0]), case[3])

    def test_level_counts_equal_class_numbers(self):
        for case in CASES:
            inst = self.instance(case)
            V = explore_from_crater(inst)
            self.assertEqual(level_counts(V), expected_level_counts(inst), case)

    def test_degrees_ell_plus_one_or_one(self):
        for case in CASES:
            inst = self.instance(case)
            V = explore_from_crater(inst)
            for j, lv in V.level.items():
                self.assertEqual(V.totals[j], case[0] + 1 if lv < case[3] else 1, (case, j, lv))

    def test_group_structure_matches_ideal_quotient(self):
        for case in CASES:
            inst = self.instance(case)
            V = explore_from_crater(inst)
            Nfac = pari.factor(inst.N)
            for j in random.Random(1).sample(sorted(V.level), min(40, len(V.level))):
                a, b = curve_for_j_fast(j, inst.p, inst.t, self.rng)
                E = pari.ellinit([a, b], inst.p)
                self.assertEqual(int(pari.ellcard(E)), inst.N)
                got = group_structure(E, inst.p, inst.N, Nfac, self.rng)
                want = predicted_group(inst.p, inst.t, case[1], case[2] * case[0] ** V.level[j])
                self.assertEqual(got, want, (case, j, V.level[j]))

    def test_group_structure_varies_with_level_when_m_positive(self):
        inst = self.instance((2, -7, 1, 6, 6, 40))
        V = explore_from_crater(inst)
        n1 = {lv: predicted_group(inst.p, inst.t, -7, 2 ** lv)[0] for lv in range(7)}
        self.assertGreater(n1[0], n1[6])
        self.assertEqual(n1[6], 1)

    def test_blind_descent_matches_bfs(self):
        for case in CASES:
            inst = self.instance(case)
            V = explore_from_crater(inst)
            k = case[3]
            for j0 in random.Random(3).sample(sorted(V.level), min(20, len(V.level))):
                r = descend_to_floor(j0, inst.p, inst.t, case[0], self.rng)
                self.assertEqual(V.level[r["floor_j"]], k, (case, j0))
                self.assertEqual(r["start_level"], V.level[j0], (case, j0))
                self.assertEqual(len(r["path"]) - 1, k - V.level[j0], (case, j0))
                # consecutive path vertices are rational ell-neighbours
                from volcano_descent import neighbors
                for u, v in zip(r["path"], r["path"][1:]):
                    self.assertIn(v, neighbors(case[0], u, inst.p)[0])

    def test_descent_from_automorphism_craters(self):
        # j = 0 (D=-3) and j = 1728 (D=-4): kernels collapse onto one j under the
        # extra automorphisms.  Regression: the crater once looked like a floor.
        for case in [(2, -3, 1, 5, 3, 40), (2, -4, 1, 5, 2, 40), (3, -4, 1, 4, 1, 40), (2, -3, 1, 5, 0, 40)]:
            inst = self.instance(case)
            V = explore_from_crater(inst)
            starts = [j for j, lv in V.level.items() if lv <= 1]
            for j0 in starts:
                r = descend_to_floor(j0, inst.p, inst.t, case[0], self.rng)
                self.assertEqual(V.level[r["floor_j"]], case[3], (case, j0))
                self.assertEqual(r["start_level"], V.level[j0], (case, j0))
                self.assertEqual(len(r["path"]) - 1, case[3] - V.level[j0], (case, j0))

    def test_local_obstruction_is_reported(self):
        # D = -7 splits 2; with f odd every element of Z + f O_K of odd norm is impossible.
        with self.assertRaises(RuntimeError):
            find_instance(11, -7, 1, 3, 0, 40, self.rng, max_tries=3000)

    def test_class_number_formula_spot(self):
        self.assertEqual(class_number(-7 * 16), 2)
        self.assertEqual(class_number(-3 * 9), 1)


class RhoTests(unittest.TestCase):
    def test_all_modes_solve(self):
        rng = random.Random(7)
        for case in [(2, -3, 1, 3, 0, 30), (2, -4, 1, 3, 0, 30), (3, -8, 1, 3, 0, 32)]:
            inst = find_instance(*case, rng, big_prime_slack=20)
            V = explore_from_crater(inst)
            Nfac = pari.factor(inst.N)
            levels = {}
            for j, lv in V.level.items():
                levels.setdefault(lv, []).append(j)
            for lv in (0, case[3]):
                j = levels[lv][0]
                a, b = curve_for_j_fast(j, inst.p, inst.t, rng)
                G, o, nf = pick_base_point(inst.p, a, b, inst.N, Nfac, rng)
                for mode in ("plain", "neg", "aut"):
                    C = Curve(inst.p, a, b)
                    x = rng.randrange(1, o)
                    Q = C.mul(x, G)
                    res = solve_ph(C, G, Q, nf, mode, rng)
                    self.assertEqual(C.mul(res["x"], G), Q, (case, lv, mode))

    def test_automorphism_orbit_is_smaller_at_crater(self):
        # sanity on the mechanism, not a measurement: mu_6 acts, so aut != neg only at j = 0
        rng = random.Random(9)
        inst = find_instance(2, -3, 1, 3, 0, 30, rng, big_prime_slack=20)
        a, b = curve_for_j_fast(0, inst.p, inst.t, rng)
        self.assertEqual(a, 0)


if __name__ == "__main__":
    unittest.main()
