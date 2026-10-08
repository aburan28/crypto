"""Transported panel replay and rank-gap accounting remain diagnostic-only."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from oracle import InvalidEvidence
from replay_paired_n17_evidence import EVIDENCE, rank_gap, replay, retained_files
from tournament import read


class PairedN17EvidenceTests(unittest.TestCase):
    def test_complete_retained_panel_replays_without_producer_execution(self):
        result = replay()
        stored = read(EVIDENCE/'REPLAY-v2.json')
        stored.pop('posthoc_replay_wall_ns')
        self.assertEqual(result, stored)
        self.assertFalse(result['promotion_eligible'])
        self.assertFalse(result['headline_online_admissible'])
        self.assertIsNone(result['online_speedup'])
        self.assertEqual(result['f5_report_status'], 'incomplete')
        self.assertEqual(result['rank_gap']['f5_rank'], 28)
        self.assertEqual([row['trial'] for row in
                          result['rank_gap']['sat_rows_resolving_gap']], [164, 173])

    def test_published_table_retains_original_keys_observations_and_failures(self):
        files = retained_files(EVIDENCE)
        panel = read(EVIDENCE/'DIAGNOSTIC-PANEL.json')
        self.assertEqual(panel['evidence_archive_sha256'],
                         read(EVIDENCE/'receipt.json')['archive_sha256'])
        self.assertEqual(len(panel['rows']), 5)
        for row in panel['rows']:
            alias = row['arm']
            original = (read(EVIDENCE/'REPLAY-v2.json')['sat_audits']['v2']
                        if alias == 'sat-v2' else
                        json.loads(files[alias+'/result.json']))
            for key in ('candidate_id', 'workload_id', 'run_id'):
                self.assertEqual(row[key], original[key])
            self.assertEqual(row['observed_online_wall_ns'],
                             original['online_wall_ns'])
            self.assertIsNone(row['admissible_online_wall_ns'])
            self.assertIsNone(row['online_speedup'])
            self.assertIsNone(row['normalized_S'])
            self.assertFalse(row['promotion_eligible'])
            if alias == 'f5':
                self.assertFalse(row['verified'])
                self.assertEqual(row['ordinary_statuses'],
                                 {'witness': 61, 'incomplete': 195})
                self.assertEqual(row['final_rank'], 28)
                self.assertEqual(row['target_attempts'], 0)

    def test_changed_archive_fails_before_any_replay(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            (path/'receipt.json').write_bytes((EVIDENCE/'receipt.json').read_bytes())
            (path/'evidence.tar.gz').write_bytes(b'changed')
            with self.assertRaisesRegex(InvalidEvidence, 'archive changed'):
                retained_files(path)

    def test_column_encoding_is_normalized_but_different_point_is_rejected(self):
        f5 = dict(columns=2, relation_matrix=dict(modulus='101',
                  column_points=[['3', '4'], ['5', '6']],
                  rows=[dict(entries=[[0, '1']], rhs='7')]),
                  collection_reports=[dict(attempts=[dict(trial=0, a=7,
                                            pdp=dict(outcome='incomplete'))])])
        sat = dict(matrix=dict(column_points=[[3, 4], [5, 6]],
                              rows=[dict(entries=[[1, '1']], scalar=7)]),
                   collection=[dict(trial=0, probe_scalar=7,
                                    status='VALID_POINT_WITNESS')])
        self.assertEqual(rank_gap(f5, sat)['nullspace_vector'], [0, 1])
        changed = copy.deepcopy(sat)
        changed['matrix']['column_points'][1] = [5, 7]
        with self.assertRaisesRegex(InvalidEvidence, 'columns differ'):
            rank_gap(f5, changed)


if __name__ == '__main__':
    unittest.main()
