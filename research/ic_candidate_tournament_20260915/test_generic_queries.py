"""Adversarial controls for query accounting, not performance admission."""
import copy
import unittest

from generic_queries import verify_queries
from oracle import Curve, InvalidEvidence
from test_certificate import FIXTURE, REPORT


def witnessed_report():
    report = copy.deepcopy(REPORT)
    report['query_schema_version'] = 1
    attempts = [dict(trial=row['trial'], a=row['a'], b=0,
                     pdp=dict(outcome='witness', points=row['points'], stats={'family': 'none'}))
                for row in report['relations']]
    report['collection_reports'] = [dict(trials=len(attempts), relations=len(attempts),
                                         attempts=attempts)]
    for solution in report['solutions']:
        rel = solution['relation']
        solution['attempts'] = [dict(trial=0, a=rel['a'], b=rel['b'],
                                     pdp=dict(outcome='witness', points=rel['points'],
                                              stats={'family': 'none'}))]
    return report


class GenericQueryTests(unittest.TestCase):
    def test_existing_certificates_replay_in_query_ledger(self):
        result = verify_queries(witnessed_report(), FIXTURE, 3)
        self.assertEqual(result['collection_queries'], 8)
        self.assertEqual(result['descent_queries'], 2)
        self.assertEqual(result['all_outcomes'], {'witness': 10})
        self.assertFalse(result['promotion_eligible'])

    def test_rejects_dropped_reordered_and_forged_queries(self):
        def attempts(report):
            return report['collection_reports'][0]['attempts']
        mutations = [
            lambda r: attempts(r).pop(),
            lambda r: attempts(r).reverse(),
            lambda r: attempts(r)[0].update(a=0),
            lambda r: attempts(r)[0]['pdp'].update(outcome='proved_unsat', points=None),
            lambda r: attempts(r)[0]['pdp'].update(outcome='identity', points=None),
            lambda r: attempts(r)[0]['pdp']['points'].__setitem__(0, 999999),
            lambda r: r['solutions'][0]['attempts'].clear(),
            lambda r: r['solutions'][0]['attempts'][0].update(a=1),
            lambda r: r['solutions'][0].update(recovered='0'),
            lambda r: r['relations'].pop(),
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                r = witnessed_report()
                mutation(r)
                with self.assertRaises(InvalidEvidence):
                    verify_queries(r, FIXTURE, 3)

    def test_failed_descent_retains_attempts_without_becoming_success(self):
        r = witnessed_report()
        s = r['solutions'][0]
        s.update(recovered=None, relation=None)
        s['attempts'][0]['pdp'].update(outcome='unresolved', points=None)
        result = verify_queries(r, FIXTURE, 3)
        self.assertEqual(result['descent_queries'], 2)
        self.assertEqual(result['all_outcomes']['unresolved'], 1)

    def test_negative_answers_are_independently_checked(self):
        curve = Curve(FIXTURE)
        r = dict(mode='ic', fixture=copy.deepcopy(FIXTURE), query_schema_version=1,
                 factor_base=[FIXTURE['generator']], trials=1, relations=[], solutions=[])
        attempt = dict(trial=0, a=1, b=0,
                       pdp=dict(outcome='proved_unsat', points=None, stats={'family': 'none'}))
        r['collection_reports'] = [dict(trials=1, relations=0, attempts=[attempt])]
        verify_queries(r, FIXTURE, 2)  # G != G+G in this nontrivial prime subgroup.
        attempt['a'] = 2
        self.assertEqual(curve.mul(curve.g, 2), curve.add(curve.g, curve.g))
        with self.assertRaisesRegex(InvalidEvidence, 'false negative'):
            verify_queries(r, FIXTURE, 2)

    def test_incomplete_counters_cannot_certify_unsat(self):
        r = witnessed_report()
        pdp = r['collection_reports'][0]['attempts'][0]['pdp']
        pdp.update(outcome='proved_unsat', points=None,
                   stats=dict(family='sat', backend='native_xor',
                              stats=dict(exhausted=True, unsupported=False, refuted=True,
                                         spurious=0, solver_calls=1)))
        with self.assertRaisesRegex(InvalidEvidence, 'incomplete refutation'):
            verify_queries(r, FIXTURE, 3)
        pdp['stats']['stats'].update(exhausted=False, spurious=1)
        with self.assertRaisesRegex(InvalidEvidence, 'invalid model'):
            verify_queries(r, FIXTURE, 3)


if __name__ == '__main__':
    unittest.main()
