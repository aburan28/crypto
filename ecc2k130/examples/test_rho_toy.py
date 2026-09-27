#!/usr/bin/env python3
import random
import unittest

import rho_toy as toy


class RhoToyTests(unittest.TestCase):
    def test_recovers_a_planted_logarithm(self):
        rng = random.Random(1)
        r, G, s = toy.setup(rng)
        planted = rng.randrange(1, r)
        Q = toy.scalar_mul(planted, G)
        k, stats, hit = toy.search(r, G, Q, s, toy.normal_masks(), rng)
        self.assertIsNotNone(k, stats)
        self.assertEqual(k, planted)
        self.assertEqual(toy.scalar_mul(k, G), Q)
        self.assertGreater(stats["dps"], 1)
        self.assertIsNotNone(hit)

    def test_main_prints_a_verified_answer(self):
        self.assertEqual(toy.main(), 0)


if __name__ == "__main__":
    unittest.main()
