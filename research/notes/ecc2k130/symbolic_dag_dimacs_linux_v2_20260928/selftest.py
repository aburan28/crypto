#!/usr/bin/env python3
"""No-network controls for the one-shot and pre-dispatch archive boundary."""
from __future__ import annotations

import contextlib
import gzip
import hashlib
import io
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import run
import v2_child
from ci_replay import (FIRST, _check_relation_stream, _replay_n131,
                       check_archive, check_static)


class GateControls(unittest.TestCase):
    def setUp(self) -> None:
        self.frozen = check_static()
        self.synthetic = {'release_pr_number': 99999,
                          'release_branch': self.frozen['release_branch']}
        self.current = {'id': 123, 'event': 'pull_request',
                        'path': '.github/workflows/' + run.WORKFLOW,
                        'workflow_id': 456,
                        'head_branch': self.synthetic['release_branch'],
                        'pull_requests': [{'number': 99999}]}

    def runs(self, rows: list[dict]) -> str:
        return '\n'.join(json.dumps(row) for row in rows) + '\n'

    def test_one_shot_accepts_current_event_only(self) -> None:
        with patch.object(run.subprocess, 'check_output',
                          side_effect=[json.dumps(self.current),
                                       self.runs([self.current])]):
            observed = run._one_shot_run(self.synthetic, 123)
        self.assertEqual(observed['matching_labeled_runs'], [123])

    def test_readded_label_or_missing_current_run_is_refused(self) -> None:
        prior = {**self.current, 'id': 122, 'pull_requests': []}
        with patch.object(run.subprocess, 'check_output',
                          side_effect=[json.dumps(self.current),
                                       self.runs([prior, self.current])]):
            with self.assertRaisesRegex(RuntimeError, 'already has a labeled run'):
                run._one_shot_run(self.synthetic, 123)
        with patch.object(run.subprocess, 'check_output',
                          side_effect=[json.dumps(self.current),
                                       self.runs([prior])]):
            with self.assertRaisesRegex(RuntimeError, 'already has a labeled run'):
                run._one_shot_run(self.synthetic, 123)

    def test_held_release_cannot_start_a_child(self) -> None:
        with patch.object(run, '_event_context', side_effect=AssertionError(
                'event context must not be inspected while held')):
            with self.assertRaisesRegex(RuntimeError, 'v2 remains held'):
                run.release_gate(self.frozen, run.git('rev-parse', 'HEAD'))

    def test_capped_child_gate_only_does_not_call_git(self) -> None:
        with patch.object(v2_child, 'git', side_effect=AssertionError(
                'child gate-only must not call Git')):
            with patch.object(sys, 'argv', ['v2_child.py', '--gate-only']):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    self.assertEqual(v2_child.main(), 0)
        self.assertEqual(json.loads(output.getvalue())['decision'],
                         'HASH_ONLY_NO_MEASURED_CHILD')

    def test_capped_child_refuses_held_dispatch(self) -> None:
        with self.assertRaisesRegex(RuntimeError, 'v2 remains held'):
            v2_child.dispatch_gate(self.frozen, run.git('rev-parse', 'HEAD'),
                                   Path('/tmp/toy'), Path('/tmp/DISPATCH.json'),
                                   '0' * 64)

    def test_capped_child_accepts_only_sealed_dispatch(self) -> None:
        released = dict(self.frozen)
        released.update(status='RELEASED',
                        release_main_head=self.frozen['base_main_head'],
                        release_pr_number=921)
        head = run.git('rev-parse', 'HEAD')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            probe = root / 'target_rlimit_probe.json'
            probe.write_text('{}\n')
            gate = {
                'reviewed_head': head, 'checkout_head': head, 'pr_head': head,
                'freeze_sha256': run.sha(run.HERE / 'FROZEN.json'),
                'release_main_head': released['release_main_head'],
                'pr_number': 921,
                'preparation_freeze_sha256':
                    released['preparation']['freeze_sha256'],
                'linux_binary_sha256': released['linux_binary_sha256'],
                'event': {'run_id': 123, 'event_head': head},
                'one_shot': {'run_id': 123, 'matching_labeled_runs': [123]},
                'target_cap_probe': {'receipt_sha256': run.sha(probe)},
            }
            dispatch = root / 'DISPATCH.json'
            dispatch.write_text(json.dumps(gate, sort_keys=True) + '\n')
            with patch.dict(v2_child.os.environ, {'GITHUB_ACTIONS': 'true',
                                                   'GITHUB_RUN_ID': '123'}):
                with patch.object(v2_child.sys, 'platform', 'linux'):
                    with patch.object(v2_child.os, 'uname',
                                      return_value=types.SimpleNamespace(
                                          machine='x86_64')):
                        with patch.object(v2_child, 'git', return_value=head):
                            self.assertEqual(
                                v2_child.dispatch_gate(released, head, root / 'toy',
                                                       dispatch, run.sha(dispatch)),
                                gate)
                            with self.assertRaisesRegex(RuntimeError,
                                                        'dispatch file path or digest'):
                                v2_child.dispatch_gate(released, head, root / 'toy',
                                                       dispatch, '0' * 64)

    def test_n131_rejects_original_equal_dimension_fake(self) -> None:
        # The former replay accepted this arbitrary one-variable CNF because
        # only its self-reported dimensions were compared with its header.
        raw = b'p cnf 1 1\n1 0\n'
        packed = gzip.compress(raw, mtime=0)
        modulus = json.loads((FIRST / 'INPUT.json').read_text())[
            'n131_single_edge_modulus']
        result = {
            'decision': 'PASS', 'single_edge_only': True, 'n': 131,
            'modulus': modulus, 'cnf_gzip_sha256': hashlib.sha256(packed).hexdigest(),
            'cnf_gzip_bytes': len(packed),
            'cnf': {'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw),
                    'variables': 1, 'clauses': 1,
                    'dag': {'variables': 920, 'model_limbs': 15, 'total_nodes': 922}},
        }
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory)
            (archive / 'n131').mkdir()
            (archive / 'n131/single_edge.cnf.gz').write_bytes(packed)
            with self.assertRaisesRegex(AssertionError, 'n131 CNF/width/node metadata'):
                _replay_n131(archive, result, self.frozen)

    def test_relation_stream_rejects_same_dimensions_mutations(self) -> None:
        sys.path.insert(0, str(FIRST))
        from export import build_relation, write_cnf
        from verify import parse_cnf

        relation = build_relation(2, 0x7)
        with tempfile.TemporaryDirectory() as directory:
            cnf = Path(directory) / 'relation.cnf'
            meta = write_cnf(relation, cnf, byte_cap=1 << 20)
            original = cnf.read_text().splitlines()
            dimensions = parse_cnf(cnf, keep_clauses=False)[:2]
            _check_relation_stream(relation, cnf, meta['clauses'])

            changed_gate = original.copy()
            literals = changed_gate[3].split()
            literals[0] = str(-int(literals[0]))
            changed_gate[3] = ' '.join(literals)
            cnf.write_text('\n'.join(changed_gate) + '\n')
            self.assertEqual(parse_cnf(cnf, keep_clauses=False)[:2], dimensions)
            with self.assertRaisesRegex(AssertionError, 'relation .* clauses changed'):
                _check_relation_stream(relation, cnf, meta['clauses'])

            changed_output = original.copy()
            changed_output[-1] = f'-{relation.output + 1} 0'
            cnf.write_text('\n'.join(changed_output) + '\n')
            self.assertEqual(parse_cnf(cnf, keep_clauses=False)[:2], dimensions)
            with self.assertRaisesRegex(AssertionError, 'output assertion changed'):
                _check_relation_stream(relation, cnf, meta['clauses'])

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
