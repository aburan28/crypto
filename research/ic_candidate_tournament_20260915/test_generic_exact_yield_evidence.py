"""Independent direct-triple control of the consequential exact-yield claim."""
import hashlib
import json
from pathlib import Path
import tarfile
import unittest

from oracle import Curve
from run_generic_exact_yield_audit import HERE, PANEL, PANEL_SHA256, digest_bytes
from tournament import read


RESULT = (HERE / 'goal_20260924/generic-exact-yield-audit/RESULT.json')
RESULT_SHA256 = 'a7119332e6b50dfd9bceb3b15dcf6d62bf579ff065bb5026e49ec0a5db09329c'


class ExactYieldEvidenceTests(unittest.TestCase):
    def test_every_exact_label_matches_direct_three_point_enumeration(self):
        panel, result = read(PANEL), read(RESULT)
        archive = HERE / panel['parent_evidence_file']
        with tarfile.open(archive, mode='r:gz') as evidence:
            for cell in panel['cells']:
                rows = [row for row in result['rows'] if row['cell'] == cell]
                reported = next(row for row in rows if 'queries' in row)
                report = json.load(evidence.extractfile(
                    f"jobs/{cell}/{reported['solver']}/stdout.json"))
                curve = Curve(report['fixture'])
                base = [curve.decode(value) for value in report['factor_base']]
                direct_sums = set()
                for i, first in enumerate(base):
                    for j in range(i, len(base)):
                        pair = curve.add(first, base[j])
                        for k in range(j, len(base)):
                            direct_sums.add(curve.add(pair, base[k]))
                for row in rows:
                    if 'queries' not in row:
                        self.assertEqual(row['analysis_status'], 'UNKNOWN_NO_REPORT')
                        continue
                    for query in row['queries']:
                        target = curve.mul(curve.g, query['a'])
                        self.assertEqual(query['exact_feasible'], target in direct_sums,
                                         (cell, row['solver'], query['trial']))
                        if query['exact_witness_indices'] is not None:
                            summands = [base[i] for i in query['exact_witness_indices']]
                            self.assertEqual(curve.add(curve.add(summands[0], summands[1]),
                                                       summands[2]), target)

    def test_n19_direct_triples_independently_confirm_sat_misses(self):
        panel = read(PANEL)
        self.assertEqual(digest_bytes(PANEL.read_bytes()), PANEL_SHA256)
        self.assertEqual(hashlib.sha256(RESULT.read_bytes()).hexdigest(), RESULT_SHA256)
        result = read(RESULT)
        self.assertEqual(result['panel_sha256'], PANEL_SHA256)
        self.assertEqual(result['parent_evidence_sha256'], panel['parent_evidence_sha256'])
        self.assertEqual(result['audited_report_rows'], 16)
        self.assertEqual(result['timeout_rows'], 4)
        self.assertEqual(len(result['rows']), 20)
        self.assertEqual(result['exact_distinct_queries'], 137)

        archive = HERE / panel['parent_evidence_file']
        self.assertEqual(hashlib.sha256(archive.read_bytes()).hexdigest(),
                         panel['parent_evidence_sha256'])
        with tarfile.open(archive, mode='r:gz') as evidence:
            report = json.load(evidence.extractfile('jobs/n19a0/f5/stdout.json'))
        curve = Curve(report['fixture'])
        base = [curve.decode(value) for value in report['factor_base']]

        # This deliberately does not call the measured PDP or the registered
        # pair-index oracle. Enumerate each three-point group sum directly.
        all_triple_sums = set()
        for i, first in enumerate(base):
            for j in range(i, len(base)):
                pair_sum = curve.add(first, base[j])
                for k in range(j, len(base)):
                    all_triple_sums.add(curve.add(pair_sum, base[k]))

        rows = {(row['cell'], row['solver']): row for row in result['rows']}
        expected = [q['a'] for q in rows['n19a0', 'f5']['queries']]
        self.assertEqual(len(expected), 24)
        exact_positive = {trial for trial, a in enumerate(expected)
                          if curve.mul(curve.g, a) in all_triple_sums}
        self.assertEqual(exact_positive, {3, 10})
        for arm in ('f5', 'sat_xor', 'sat_cnf'):
            row = rows['n19a0', arm]
            self.assertEqual([q['a'] for q in row['queries']], expected)
            self.assertEqual({q['trial'] for q in row['queries'] if q['exact_feasible']},
                             exact_positive)
            self.assertEqual(row['feasible_but_incomplete'], 0 if arm == 'f5' else 2)
        self.assertEqual([rows['n19a0', arm]['observed_witnesses']
                          for arm in ('f5', 'sat_xor', 'sat_cnf')], [2, 0, 0])


if __name__ == '__main__':
    unittest.main()
