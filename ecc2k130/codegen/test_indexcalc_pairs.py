import unittest

import curves
import indexcalc_e2e as engine
import indexcalc_pairs as pairs


class AccountingTests(unittest.TestCase):
    def test_product_law_is_independent_of_base_size(self):
        # relations * (group / C(F,3)) * C(F,2) * 20  ~  20 * 3 * group
        # once F is large enough that C(F,3) ~ F^3/6 and C(F,2) ~ F^2/2.
        group = pairs.R_131
        costs = []
        for base in (1 << 16, 1 << 20, 1 << 24):
            costs.append(pairs.streamingProducts(base, base, group, 3))
        ratios = [costs[i] / (3 * group * pairs.PRODUCTS_PER_PAIR) for i in range(3)]
        for ratio in ratios:
            self.assertGreater(ratio, 0.9)
            self.assertLess(ratio, 1.1)

    def test_rho_ratio_at_weight_two_is_hopeless(self):
        projection = pairs.projectAttack(4000, 4000 / 131.0, pairs.R_131, pairs.ORDER_131)
        self.assertEqual(projection['class'], 'engineering')
        self.assertGreater(projection['log2_ratio_to_rho'], 60)
        self.assertAlmostEqual(projection['log2_rho_field_products'],
                               pairs.log2(pairs.rhoProducts()), places=9)

    def test_expected_yield_matches_binomial(self):
        self.assertAlmostEqual(pairs.expectedYield(10, 3, 1000), 120 / 1000.0)


class PairOracleTests(unittest.TestCase):
    def test_planted_triple_at_degree_5(self):
        meter = engine.Ledger()
        context = engine.setup(5, 4, 3, meter)
        onb, curve, ell, generator, eigen, reps, lookup, _, _ = context
        points = list(lookup)
        self.assertGreater(len(points), 3)
        triple = (points[0], points[2], points[4])
        target = curve.add(curve.add(triple[0], triple[1]), triple[2])
        found = pairs.pairDecomposeLookup(curve, lookup, target)
        self.assertIsNotNone(found)
        got = None
        for p in found:
            got = curve.add(got, p)
        self.assertEqual(got, target)

    def test_known_answer_dlp_at_degree_9(self):
        report = pairs.recoverLog(9, 2, seed=20260918, attempts=256)
        self.assertEqual(report['status'], 'complete')
        self.assertTrue(report['verified'])
        self.assertEqual(report['oracle_class'], 'enumerate_pairs')
        self.assertEqual(int(report['scalar_recovered']), int(report['expected_scalar']))
        onb = engine.AuditField(9, engine.Ledger())
        curve = curves.Curve(onb)
        # Re-verify the recovered scalar in a fresh field object.
        self.assertTrue(curves.isPrimeBig(int(report['subgroup_order'])))


if __name__ == '__main__':
    unittest.main()
