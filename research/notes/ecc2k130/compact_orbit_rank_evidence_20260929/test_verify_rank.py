"""Mutation controls for independent rank replay."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from verify_rank import verify

HERE = Path(__file__).resolve().parent
RAW = HERE / 'evidence/local_control_20260929/traced'


class RankReplayMutationTests(unittest.TestCase):
    def test_archived_control_passes(self):
        result = verify(RAW / 'rank.jsonl', RAW / 'base.jsonl', RAW / 'summary.jsonl')
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['rank'], 2)
        self.assertEqual(result['representative_logs_verified'], 2)

    def reject_mutation(self, mutate):
        rows = [json.loads(line) for line in (RAW / 'rank.jsonl').read_text().splitlines()]
        mutate(rows)
        with tempfile.TemporaryDirectory() as folder:
            trace = Path(folder) / 'mutated.jsonl'
            trace.write_text(''.join(json.dumps(row, sort_keys=True) + '\n' for row in rows))
            with self.assertRaises(AssertionError):
                verify(trace, RAW / 'base.jsonl', RAW / 'summary.jsonl')

    def test_wrong_group_relation_rejected(self):
        self.reject_mutation(lambda rows: rows[1]['target'].__setitem__(0, rows[1]['target'][0] ^ 1))

    def test_wrong_matrix_row_rejected(self):
        self.reject_mutation(lambda rows: rows[1]['row'].__setitem__(0, rows[1]['row'][0] ^ 1))

    def test_wrong_base_log_rejected(self):
        self.reject_mutation(lambda rows: rows[-1]['logs'].__setitem__(0, rows[-1]['logs'][0] ^ 1))


if __name__ == '__main__':
    unittest.main()
