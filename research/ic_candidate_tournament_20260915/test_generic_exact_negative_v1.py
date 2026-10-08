"""Adversarial exact-negative controls, including geometric torsion points."""
import itertools
import unittest

from generic_exact_negative_v1 import ExactThreeSum
from oracle import Curve, InvalidEvidence
from test_certificate import FIXTURE


class ExactNegativeTests(unittest.TestCase):
    def setUp(self):
        self.curve = Curve(FIXTURE)

    def test_agrees_with_direct_triple_enumeration_and_keeps_repeated_indices(self):
        c = self.curve
        for base in ([c.g], [c.g, c.neg(c.g)], [(0, 1), c.g],
                     [c.g, c.mul(c.g, 2), c.mul(c.g, 7)]):
            with self.subTest(base=base):
                proof = ExactThreeSum(c, base, 3)
                triples = {c.add(c.add(p, q), r)
                           for p, q, r in itertools.product(base, repeat=3)}
                for target in triples - {None}:
                    with self.assertRaisesRegex(InvalidEvidence, 'false negative'):
                        proof.verify_absent(target, stage='collection', trial=0, a=1, b=0)
                for scalar in range(1, 40):
                    target = c.mul(c.g, scalar)
                    if target not in triples:
                        proof.verify_absent(target, stage='target', trial=scalar, a=scalar, b=1)
                self.assertGreater(len(proof.receipt()['rows']), 0)
                self.assertEqual(proof.receipt()['repetition_policy'], 'all-indices-with-replacement')

    def test_torsion_and_identity_pair_sums_cannot_be_filtered_from_negative_proof(self):
        c = self.curve
        torsion = (0, 1)
        self.assertIsNotNone(c.mul(torsion, c.r))
        self.assertIsNone(c.add(torsion, torsion))
        proof = ExactThreeSum(c, [torsion, c.g], 3)
        self.assertTrue(proof.receipt()['identity_pair_sums_included'])
        # T+T+G=G, although T vanishes under cofactor projection.
        with self.assertRaisesRegex(InvalidEvidence, 'false negative'):
            proof.verify_absent(c.g, stage='collection', trial=0, a=1, b=0)

    def test_oversized_wrong_summands_duplicate_and_identity_bases_fail_closed(self):
        c = self.curve
        for base, summands in [([c.g], 2), ([c.g], True), ([c.g]*129, 3),
                               ([c.g, c.g], 3), ([None], 3), ([], 3)]:
            with self.subTest(base_size=len(base), summands=summands):
                with self.assertRaises(InvalidEvidence):
                    ExactThreeSum(c, base, summands)
        c.n = 19
        with self.assertRaisesRegex(InvalidEvidence, 'domain'):
            ExactThreeSum(c, [c.g], 3)


if __name__ == '__main__':
    unittest.main()
