#!/usr/bin/env python3
"""No-network controls for the one-shot and pre-dispatch archive boundary."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run
from ci_replay import check_archive, check_static


class GateControls(unittest.TestCase):
    def setUp(self) -> None:
        self.frozen = check_static()
        self.synthetic = {'release_pr_number': 99999,
                          'release_branch': self.frozen['release_branch']}
        self.current = {'id': 123, 'event': 'pull_request',
                        'head_branch': self.synthetic['release_branch'],
                        'pull_requests': [{'number': 99999}]}

    def runs(self, rows: list[dict]) -> str:
        return '\n'.join(json.dumps(row) for row in rows) + '\n'

    def test_one_shot_accepts_current_event_only(self) -> None:
        with patch.object(run.subprocess, 'check_output',
                          return_value=self.runs([self.current])):
            observed = run._one_shot_run(self.synthetic, 123)
        self.assertEqual(observed['matching_labeled_runs'], [123])

    def test_readded_label_or_missing_current_run_is_refused(self) -> None:
        prior = {**self.current, 'id': 122, 'pull_requests': []}
        with patch.object(run.subprocess, 'check_output',
                          return_value=self.runs([prior, self.current])):
            with self.assertRaisesRegex(RuntimeError, 'already has a labeled run'):
                run._one_shot_run(self.synthetic, 123)
        with patch.object(run.subprocess, 'check_output',
                          return_value=self.runs([prior])):
            with self.assertRaisesRegex(RuntimeError, 'already has a labeled run'):
                run._one_shot_run(self.synthetic, 123)

    def test_held_release_cannot_start_a_child(self) -> None:
        with patch.object(run, '_event_context', side_effect=AssertionError(
                'event context must not be inspected while held')):
            with self.assertRaisesRegex(RuntimeError, 'v2 remains held'):
                run.release_gate(self.frozen, run.git('rev-parse', 'HEAD'))

    def test_pre_dispatch_refusal_is_replayable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            reason = {'decision': 'PRE_DISPATCH_REFUSAL',
                      'reason': 'HELD', 'attempts_started': 0}
            (out / 'PRE_DISPATCH_REFUSAL.json').write_text(
                json.dumps(reason, sort_keys=True) + '\n')
            receipt = run._manifest_and_receipt(
                out, self.frozen, None, [], 'NOT_ADMITTED: HELD', 0.01)
            self.assertEqual(receipt['decision'], 'FAIL_OR_CENSORED')
            self.assertEqual(check_archive(out / 'receipt.json', self.frozen),
                             {'decision': 'ARCHIVED_PRE_DISPATCH_REFUSAL',
                              'phases': 0})


if __name__ == '__main__':
    unittest.main()
