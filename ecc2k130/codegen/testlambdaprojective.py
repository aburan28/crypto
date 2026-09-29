"""lambda-projective formulas, invariance and pricing (LAMBDA-PROJECTIVE.md)."""
import os
import random
import unittest

import curves
import field
import lambda_projective as lp


class LambdaProjectiveTests(unittest.TestCase):
    def test_formulas_agree_with_affine_at_published_counts(self):
        rng = random.Random(7)
        mixed, full, affine, ld = lp.checkFormulas(23, 40, rng)
        self.assertEqual(mixed, (8, 2, 0))
        self.assertEqual(full, (11, 2, 0))
        self.assertEqual(affine, (6, 2, 1))
        self.assertEqual(ld, (8, 5, 0))

    def test_twist_identity(self):
        # f(c R) is a representative of the point f(R) whenever the selector
        # is scale-invariant: the walk descends to points.
        rng = random.Random(11)
        onb = field.Onb(23)
        curve = curves.Curve(onb)
        sel = lp.invariantSelector(onb)
        for _ in range(10):
            p = lp.randomAffinePoint(curve, rng)
            r = lp.toLambdaProjective(onb, lp.toLambdaAffine(onb, p), onb.one())
            c = onb.randomElement(rng) or onb.one()
            s = (onb.mul(r[0], c), onb.mul(r[1], c), onb.mul(r[2], c))
            for _ in range(4):
                r = lp.walkStep(onb, r, sel)
                s = lp.walkStep(onb, s, sel)
                self.assertTrue(lp.sameProjectivePoint(onb, r, s))

    def test_representative_selector_splits_trails(self):
        rng = random.Random(13)
        s = lp.splitStatistics(23, 60, 8, rng)
        self.assertEqual(s['invariantSplits'], 0)
        self.assertGreater(s['splitWithin2'], 0.8)
        self.assertGreater(s['dpSpurious'], 0.0)

    def test_price_table_orders_the_variants(self):
        rows = dict((name, (alu, cl)) for name, mu, sq, iv, alu, cl in lp.priceTable((8, 2, 0), (11, 2, 0)))
        affine16 = rows['affine, Montgomery batch 16']
        projective = rows['lambda-projective mixed (table walk), hash free']
        # the affine floor of ITERATION-FUNCTION.md, re-priced with the
        # THROUGHPUT-20B routine costs, lands where that note put it
        self.assertAlmostEqual(affine16[1], 38.1, places=1)
        self.assertTrue(1090 <= affine16[0] <= 1180)
        # and lambda-projective with a free invariant hash is above it on both pipes
        self.assertGreater(projective[1], affine16[1])
        self.assertGreater(projective[0], affine16[0])
        # on the carry-less column affine is ahead from batch 2; on the
        # binding-pipe ceiling the crossover sits between batch 2 and 3
        self.assertLess(rows['affine, Montgomery batch 2'][1], projective[1])
        ceiling = lambda row: lp.ceilingBs(row[1], row[0])[1]
        self.assertLess(ceiling(rows['affine, Montgomery batch 2']), ceiling(projective))
        self.assertGreater(ceiling(rows['affine, Montgomery batch 3']), ceiling(projective))
        # and once the hash is paid by an inversion, batch-1 affine beats it on both pipes
        paid = rows['lambda-projective mixed + invariant hash by one inversion']
        self.assertLess(rows['affine, Montgomery batch 1'][0], paid[0])
        self.assertLess(rows['affine, Montgomery batch 1'][1], paid[1])

    def test_price_table_matches_the_note(self):
        # LAMBDA-PROJECTIVE.md carries the table this script prints, and the
        # campaign status page is pinned to the note; so a price change here
        # (a routine cost, a formula count) must reach the note in the same
        # commit or this fails. Every row's carry-less and ALU figure, as the
        # script formats it, has to appear as a table cell of the note.
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, '..', 'LAMBDA-PROJECTIVE.md'), encoding='utf-8') as fh:
            note = fh.read().replace('*', '')
        for name, mu, sq, iv, alu, cl in lp.priceTable((8, 2, 0), (11, 2, 0)):
            self.assertIn('| %.1f |' % cl, note, '%s: %.1f clmad is not in the note' % (name, cl))
            self.assertIn('| {:,.0f} |'.format(alu), note, '%s: %.0f ALU is not in the note' % (name, alu))


if __name__ == '__main__':
    unittest.main()
