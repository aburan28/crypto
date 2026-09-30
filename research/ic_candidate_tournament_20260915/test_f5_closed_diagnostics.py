"""Disclosed correctness controls, never natural-yield or solver measurements."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from f5_closed_diagnostics import analyze, chain_control, counters, FLAGS, INTEGER_COUNTERS
from oracle import Curve, InvalidEvidence

FIXTURE = dict(degree=17, curve_a=1, subgroup_order='65587', group_order='131174',
               cofactor='2', generator=['43693', '23339'], **{'lambda': '17184'},
               irreducible=dict(degree=17, low_terms=[0, 3]))


class F5ClosedDiagnosticsTests(unittest.TestCase):
    def test_repeated_sign_and_identity_intermediate_controls(self):
        curve = Curve(FIXTURE)
        base = [curve.g, curve.neg(curve.g), curve.mul(curve.g, 3), (0, 1)]
        controls = [([0, 0, 2], curve.mul(curve.g, 5)),
                    ([0, 1, 2], curve.mul(curve.g, 3)),
                    ([3, 3, 0], curve.g)]
        for indices, target in controls:
            result = chain_control(curve, base, indices, target)
            self.assertTrue(result['full_point_replay'])
            self.assertTrue(result['finite_chain_exists'])
            self.assertEqual(len(result['orderings']), 6)
            self.assertFalse(result['implemented_boolean_anf_validated'])
        inverse_pair = chain_control(curve, base, [0, 1, 2], base[2])
        identities = [row for row in inverse_pair['orderings'] if row['status'] == 'IDENTITY_INTERMEDIATE']
        self.assertEqual(len(identities), 2)
        self.assertTrue(all(row['s3_values'] is None and row['intermediate'] is None
                            and not row['represented_by_finite_chain'] for row in identities))

    def test_two_torsion_degeneracy_is_preserved_and_bad_witness_rejected(self):
        curve = Curve(FIXTURE)
        base = [curve.g, (0, 1)]
        result = chain_control(curve, base, [1, 1, 1], base[1])
        self.assertFalse(result['finite_chain_exists'])
        self.assertTrue(all(row['status'] == 'IDENTITY_INTERMEDIATE' for row in result['orderings']))
        with self.assertRaisesRegex(InvalidEvidence, 're-add'):
            chain_control(curve, base, [0, 0, 0], curve.mul(curve.g, 2))
        with self.assertRaisesRegex(InvalidEvidence, 'malformed'):
            chain_control(curve, base, [0, 0, True], curve.mul(curve.g, 3))

    def attempt(self, outcome='incomplete'):
        stats = {key: 0 for key in INTEGER_COUNTERS}
        stats.update(exhausted=True, unsupported=False, reductions=4096, splits=1573,
                     max_degree_built=3, propagations=6569)
        return dict(trial=0, a=8015, b=0, pdp=dict(outcome=outcome, points=None,
                    stats=dict(family='groebner', engine={'MatrixF5': {'max_degree': 3}}, stats=stats)))

    def test_reduction_budget_and_split_effort_remain_distinct(self):
        result = counters([self.attempt()], dict(node_budget=4096, groebner_degree=3),
                          sorted(INTEGER_COUNTERS | FLAGS))
        row = result['ledger'][0]
        self.assertTrue(row['at_reduction_budget'])
        self.assertEqual(row['stats']['reductions'], 4096)
        self.assertEqual(row['stats']['splits'], 1573)
        self.assertEqual(result['by_outcome']['incomplete']['attempts'], 1)
        self.assertIn('cache_hits_and_misses', result['unresolved_instrumentation'])
        with self.assertRaisesRegex(InvalidEvidence, 'completion flags'):
            counters([self.attempt('proved_unsat')], dict(node_budget=4096, groebner_degree=3),
                     sorted(INTEGER_COUNTERS | FLAGS))

    def test_missing_counter_boolean_counter_or_changed_dispatch_is_rejected(self):
        config = dict(node_budget=4096, groebner_degree=3)
        fields = sorted(INTEGER_COUNTERS | FLAGS)
        for mutate in (
            lambda row: row['pdp']['stats']['stats'].pop('reductions'),
            lambda row: row['pdp']['stats']['stats'].update(reductions=True),
            lambda row: row['pdp']['stats'].update(engine={'MatrixF4': {'max_degree': 3}}),
        ):
            row = copy.deepcopy(self.attempt())
            mutate(row)
            with self.assertRaises(InvalidEvidence):
                counters([row], config, fields)

    def test_wrong_external_archive_is_rejected_before_reading_members(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'receipt.json').write_text(json.dumps(dict(archive_sha256='0'*64, archive_bytes=1)))
            with self.assertRaisesRegex(InvalidEvidence, 'externally fixed'):
                analyze(root)


if __name__ == '__main__':
    unittest.main()
