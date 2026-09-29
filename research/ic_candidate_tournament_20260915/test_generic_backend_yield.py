"""Regression controls for censored natural-query and repeated-point accounting."""
import unittest

from generic_backend_yield import summarize
from oracle import InvalidEvidence


def row(case, repetition, *, witness, queries, status='VERIFIED'):
    audited = queries is not None
    return dict(stage='development', arm='generic_f4_dense', cell='n17a1',
        case=case, repetition=repetition, execution_status=status,
        audited=audited, query_count=queries, witness_count=witness if audited else None,
        outcome_mix={'witness':witness, 'unresolved':queries-witness} if audited else None,
        accepted_rows=witness if audited else None, rank=1 if witness else 0,
        base_points=272 if audited else None, folded_columns=8 if audited else None,
        censor_reason=None if audited else 'no complete profiler report',
        complete_instruction_phases=None, target_descent_queries=0 if audited else None,
        target_outcome_mix={} if audited else None)


class NaturalYieldTests(unittest.TestCase):
    def test_repetitions_do_not_become_new_points_and_censoring_is_explicit(self):
        observations = [row('p0', rep, witness=1, queries=2) for rep in range(3)]
        observations += [row('p1', rep, witness=0, queries=None,
                             status='TIMEOUT') for rep in range(3)]
        result = summarize(observations)[0]
        self.assertEqual(result['observed_ordinary_queries'], 6)
        self.assertEqual(result['verified_witness_queries'], 3)
        self.assertEqual(result['censored_runs'], 3)
        self.assertEqual(result['natural_witness_rate']['distinct_points'], 1)
        self.assertEqual(result['natural_witness_rate']['rate'], .5)
        self.assertEqual(result['qualification'], 'incomplete')

    def test_changed_query_history_across_process_repetitions_is_rejected(self):
        observations = [row('p0', rep, witness=1, queries=2) for rep in range(3)]
        observations[2]['outcome_mix'] = {'witness':1, 'unresolved':0, 'identity':1}
        with self.assertRaises(InvalidEvidence):
            summarize(observations)

    def test_zero_yield_retains_nonzero_distinct_point_uncertainty(self):
        observations = [row(f'p{point}', rep, witness=0, queries=8)
                        for point in range(3) for rep in range(3)]
        rate = summarize(observations)[0]['natural_witness_rate']
        self.assertEqual(rate['rate'], 0)
        self.assertEqual(rate['ci95'], [0, 0])
        self.assertGreater(rate['point_mean_hoeffding95'][1], 0)

    def test_single_process_schedule_uses_distinct_points_only(self):
        observations = [row(f'p{i}', 0, witness=i % 2, queries=2)
                        for i in range(4)]
        result = summarize(observations, repetitions=1)[0]
        self.assertEqual(result['natural_witness_rate']['distinct_points'], 4)
        self.assertEqual(result['natural_witness_rate']['rate'], .25)
        self.assertEqual(result['scheduled_runs'], 4)
        with self.assertRaises(InvalidEvidence):
            summarize(observations)

    def test_audited_zero_query_point_is_not_zero_natural_yield(self):
        observations = [row('p0', 0, witness=0, queries=0)]
        rate = summarize(observations, repetitions=1)[0]['natural_witness_rate']
        self.assertIsNone(rate['rate'])
        self.assertEqual(rate['distinct_points'], 0)
        self.assertEqual(rate['zero_query_points'], 1)

    def test_zero_query_point_prevents_full_point_law_uncertainty_claim(self):
        observations = [row('p0', 0, witness=0, queries=0),
                        row('p1', 0, witness=1, queries=2)]
        rate = summarize(observations, repetitions=1)[0]['natural_witness_rate']
        self.assertEqual(rate['rate'], .5)
        self.assertEqual(rate['zero_query_points'], 1)
        self.assertIsNone(rate['point_mean_hoeffding95'])


if __name__ == '__main__':
    unittest.main()
