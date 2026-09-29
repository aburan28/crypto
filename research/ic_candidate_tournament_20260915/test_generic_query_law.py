"""Independent RNG vectors and valid-group/invalid-sampler adversarial controls."""
import copy
import itertools
import json
from pathlib import Path
import unittest

from generic_queries import verify_queries
from generic_query_law import StdRng08, descent_coefficients, verify_query_law, verify_rust_vectors
from oracle import InvalidEvidence

ROOT = Path(__file__).parent
ARCHIVE = ROOT / 'goal_20260924/generic-public-inputs/final-worker-raw.jsonl'


def archived():
    return [row for line in ARCHIVE.read_text().splitlines()
            if (row := json.loads(line))['job']['mode'] == 'ic']


class GenericQueryLawTests(unittest.TestCase):
    def test_frozen_fresh_controls_and_window_dispatch(self):
        path = ROOT / 'goal_20260924/generic-query-law/worker-raw.jsonl'
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        self.assertEqual(len(rows), 47)
        for row in rows:
            with self.subTest(case=row['name']):
                verify_query_law(row['report'], row['report']['fixture'], row['job'])
        row = next(row for row in rows if row['name'] == 'n9-window1-dense')
        row['job']['config']['collection_window'] = 0
        verify_queries(row['report'], row['report']['fixture'], 3)
        with self.assertRaisesRegex(InvalidEvidence, 'collection coefficient mismatch'):
            verify_query_law(row['report'], row['report']['fixture'], row['job'])

    def test_upstream_stdrng_construction_vector(self):
        # rand 0.8.8 src/rngs/std.rs::test_stdrng_construction.
        rng = StdRng08([1, 23, 456, 7890, 0, 0, 0, 0])
        self.assertEqual(rng.next_u64(), 10719222850664546238)
        child = StdRng08([rng.next_u32() for _ in range(8)])
        self.assertEqual(child.next_u64(), 14064965282130556830)

    def test_pinned_rust_streams_rejections_and_wrapping_probes(self):
        vectors = json.loads((ROOT / 'goal_20260924/generic-query-law/rust-vectors.json').read_text())
        self.assertEqual(verify_rust_vectors(vectors)['values_verified'], 7436)
        vectors['streams'][0]['values'][128] ^= 1
        with self.assertRaisesRegex(InvalidEvidence, 'differs from Rust'):
            verify_rust_vectors(vectors)

    def test_all_archived_queries_match_declared_jobs(self):
        rows = archived()
        self.assertEqual(len(rows), 35)
        for row in rows:
            with self.subTest(case=row['name']):
                receipt = verify_query_law(row['report'], row['report']['fixture'], row['job'])
                self.assertFalse(receipt['promotion_eligible'])

    def test_wrong_seed_is_rejected_without_changing_valid_group_evidence(self):
        row = archived()[0]
        verify_queries(row['report'], row['report']['fixture'], 2)
        row['job']['algorithm_seed'] ^= 1
        with self.assertRaisesRegex(InvalidEvidence, 'collection coefficient mismatch'):
            verify_query_law(row['report'], row['report']['fixture'], row['job'])

    def test_true_replacement_witness_cannot_change_the_collection_sampler(self):
        row = archived()[0]
        report = row['report']
        attempts = report['collection_reports'][0]['attempts']
        attempts[0] = dict(copy.deepcopy(attempts[1]), trial=0)
        report['relations'][0] = dict(trial=0, a=attempts[0]['a'], points=attempts[0]['pdp']['points'])
        verify_queries(report, report['fixture'], 2)  # The replacement group equation is true.
        with self.assertRaisesRegex(InvalidEvidence, 'collection coefficient mismatch at 0'):
            verify_query_law(report, report['fixture'], row['job'])

    def test_true_replacement_descent_cannot_change_the_sampler(self):
        row = archived()[0]
        report = row['report']
        solution = report['solutions'][0]
        last = solution['attempts'][-1]
        order, log = int(report['fixture']['subgroup_order']), int(solution['recovered'])
        b = last['b'] % (order - 1) + 1
        a = (last['a'] - log * (b - last['b'])) % order
        last.update(a=a, b=b)
        solution['relation'].update(a=a, b=b)
        verify_queries(report, report['fixture'], 2)  # Same query point and recovered scalar.
        with self.assertRaisesRegex(InvalidEvidence, 'descent coefficient mismatch'):
            verify_query_law(report, report['fixture'], row['job'])

    def test_rejects_wrong_batch_limits_fields_targets_and_method(self):
        mutations = (
            lambda j: j['config'].update(batch_trials=4),
            lambda j: j['config'].update(max_trials=4),
            lambda j: j['config'].update(summands=3),
            lambda j: j.update(algorithm_seed=True),
            lambda j: j.update(degree=13),
            lambda j: j.update(public_targets=[]),
            lambda j: j['config'].update(solver='unknown'),
        )
        for mutate in mutations:
            row = archived()[0]
            mutate(row['job'])
            with self.assertRaises(InvalidEvidence):
                verify_query_law(row['report'], row['report']['fixture'], row['job'])

    def test_two_summand_window_is_inactive_and_walk_keeps_identities(self):
        row = archived()[0]
        row['job']['config']['collection_window'] = 1
        receipt = verify_query_law(row['report'], row['report']['fixture'], row['job'])
        self.assertEqual(receipt['collection_law'], 'sampled')
        walked = list(itertools.islice(descent_coefficients(2026092556, 31, True), 129))
        self.assertTrue(any(a == 0 for a, _ in walked))
        self.assertEqual(len({b for _, b in walked}), 1)
        for t in range(65):
            self.assertEqual(walked[t + 64], ((walked[t][0] + 1) % 31, walked[t][1]))

    def test_preparation_failure_cannot_skip_paid_trials(self):
        row = next(row for row in archived() if row['name'] == 'incomplete-pair_table')
        row['job']['config']['max_trials'] = 2
        verify_queries(row['report'], row['report']['fixture'], 2)
        with self.assertRaisesRegex(InvalidEvidence, 'failed preparation stopped early'):
            verify_query_law(row['report'], row['report']['fixture'], row['job'])


if __name__ == '__main__':
    unittest.main()
