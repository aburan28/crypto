import copy
import itertools
import json
import importlib.util
import io
from pathlib import Path
import subprocess
import tempfile
import tarfile
import unittest
from unittest.mock import patch

from f5_boolean_control import (add, audit, checked_rows, equations, evaluate, multiply,
                                products, rename_mask, rename_rows, row_space, s3_polynomial)
from oracle import Curve, InvalidEvidence


class TinyField:
    n = 3

    @staticmethod
    def fm(a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & 8:
                a ^= 11
        return out


def field_evaluate(poly, assignment):
    out = 0
    for mask, coefficient in poly.items():
        if mask & assignment == mask:
            out ^= coefficient
    return out


class BooleanControlTests(unittest.TestCase):
    def test_symbolic_s3_agrees_on_entire_small_cube(self):
        c = TinyField()
        xs = [{1 << (2*i): 1, 1 << (2*i+1): 2} for i in range(3)]
        polynomial = s3_polynomial(c, *xs)
        for assignment in range(64):
            x, y, z = [field_evaluate(poly, assignment) for poly in xs]
            xy = c.fm(x, y)
            value = c.fm(c.fm(x ^ y, x ^ y), c.fm(z, z)) ^ c.fm(xy, z) ^ c.fm(xy, xy) ^ 1
            self.assertEqual(field_evaluate(polynomial, assignment), value)
        self.assertEqual(multiply(c, {1: 1, 2: 1}, {1: 1, 2: 1}), {1: 1, 2: 1})
        self.assertEqual(add({1: 7}, {1: 7}), {})

    def test_chain_coordinates_and_permutation_preserve_all_assignments(self):
        c = TinyField()
        rows = equations(c, [1], 5)
        permutation = [2, 0, 4, 1, 5, 3]
        renamed = rename_rows(rows, permutation)
        for assignment in range(64):
            self.assertEqual([evaluate(row, assignment) for row in rows],
                             [evaluate(row, rename_mask(assignment, permutation)) for row in renamed])
        self.assertTrue(any(evaluate(row, 0) for row in rows))

    def test_row_space_is_canonical_and_faults_change_it(self):
        rows = [[1, 2], [2, 4], [1, 4]]
        columns = [1, 2, 4]
        expected = row_space(rows, columns)
        for ordering in itertools.permutations(rows):
            self.assertEqual(row_space(ordering, columns), expected)
        self.assertNotEqual(row_space([[1, 2], [2, 4], [1]], columns), expected)
        with self.assertRaisesRegex(InvalidEvidence, 'outside'):
            row_space([[8]], columns)

    def test_macaulay_products_vanish_on_roots_and_cancel_collisions(self):
        original = [[1, 2], [4, 0]]
        shifted = products(original, 3, 2)
        self.assertEqual(products([[1, 3]], 2, 3).count([1, 3]), 2)
        roots = [a for a in range(8) if all(evaluate(row, a) == 0 for row in original)]
        self.assertEqual(len(roots), 2)
        self.assertTrue(all(evaluate(row, a) == 0 for row in shifted for a in roots))

    def test_malformed_boolean_masks_rejected(self):
        for rows in ([[True]], [[1, 1]], [[8]], [[-1]]):
            with self.assertRaises(InvalidEvidence):
                checked_rows(rows, 3)

    def test_complete_control_rejects_encoding_permutation_matrix_and_witness_faults(self):
        inputs = json.loads((Path(__file__).parent/'goal_20260924/f5-boolean-system-control-20260930/inputs.json').read_text())
        curve = Curve(inputs['fixture'])
        permutation = list(range(12))+list(range(29, 35))+list(range(12, 29))
        controls = []
        for item in inputs['controls']:
            rows = equations(curve, inputs['basis'], item['target'][0])
            renamed = rename_rows(rows, permutation)
            controls.append(dict(trial=item['trial'], target=item['target'], n_vars=35, ell=6, m=3,
                direct=rows, template=rows, reused=rows, repeated=rows,
                permutation=permutation, interleaved=renamed,
                matrices=[dict(layout=layout, engine=engine, status='REDUCED',
                    rows=products(system, 35, 3), elimination_and_criterion_word_xors=0)
                    for layout, system in [('original', rows), ('interleaved', renamed)]
                    for engine in ('f4', 'f5')]))
        exported = dict(schema_version=1, basis=inputs['basis'], controls=controls)
        self.assertEqual(audit(inputs, exported)['status'], 'AUDITED_DISCLOSED_BOOLEAN_SYSTEM_CONTROL')
        changed = copy.deepcopy(exported)
        changed['controls'][0]['template'] = copy.deepcopy(changed['controls'][0]['template'])
        changed['controls'][0]['template'][0] = sorted(set(changed['controls'][0]['template'][0]) ^ {0})
        with self.assertRaisesRegex(InvalidEvidence, 'coefficient identity'):
            audit(inputs, changed)
        changed = copy.deepcopy(exported)
        changed['controls'][0]['permutation'][0:2] = [1, 0]
        with self.assertRaisesRegex(InvalidEvidence, 'permutation'):
            audit(inputs, changed)
        changed = copy.deepcopy(exported)
        changed['controls'][0]['matrices'][1]['rows'] = []
        with self.assertRaisesRegex(InvalidEvidence, 'row space differs'):
            audit(inputs, changed)
        changed_inputs = copy.deepcopy(inputs)
        changed_inputs['controls'][0]['points'][0][1] ^= 1
        with self.assertRaises(InvalidEvidence):
            audit(changed_inputs, exported)

    def test_optional_tool_version_failure_is_retained_and_first_registration_is_closed(self):
        directory = Path(__file__).parent/'goal_20260924/f5-boolean-system-control-20260930'
        spec = importlib.util.spec_from_file_location('f5_boolean_control_runner', directory/'run_control.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        version = subprocess.CompletedProcess(['ar', '--version'], 1, '', 'usage: ar')
        with patch.object(module.subprocess, 'run', return_value=version):
            receipt = module.tool_receipt(str(directory/'export.rs'), directory, {})
        self.assertEqual(receipt['exit_code'], 1)
        self.assertEqual(receipt['stderr'], 'usage: ar')
        self.assertFalse(receipt['version_probe_succeeded'])
        with self.assertRaisesRegex(InvalidEvidence, 'first registration closed'):
            module.run(directory, directory/'unused-first-registration-output', directory/'protocol.json')
        with self.assertRaisesRegex(InvalidEvidence, 'native control registration closed'):
            module.run(directory, directory/'unused-second-registration-output', directory/'protocol-v2.json')

    def test_closed_native_transport_replays_and_rejects_receipt_and_archive_faults(self):
        from replay_f5_boolean_control import replay, replay_files
        from replay_paired_n17_evidence import retained_files
        bundle = Path(__file__).parent/'goal_20260924/f5-boolean-system-control-20260930/native-v2'
        files = retained_files(bundle)
        self.assertEqual(replay_files(files)['status'], 'AUDITED_DISCLOSED_BOOLEAN_SYSTEM_CONTROL')
        for role in ('diagnostic', 'native.stdout'):
            changed = dict(files)
            changed[role] += b'changed'
            with self.assertRaises(InvalidEvidence):
                replay_files(changed)
        changed = dict(files)
        data = io.BytesIO()
        with tarfile.open(fileobj=data, mode='w:gz') as tar:
            info = tarfile.TarInfo('Cargo.toml')
            info.size = 7
            tar.addfile(info, io.BytesIO(b'changed'))
        changed['root-source.tar.gz'] = data.getvalue()
        with self.assertRaisesRegex(InvalidEvidence, 'compiled manifest'):
            replay_files(changed)
        changed = dict(files)
        process = json.loads(changed['native-process.json'])
        process['argv'][-1] += '-changed'
        changed['native-process.json'] = json.dumps(process).encode()
        with self.assertRaisesRegex(InvalidEvidence, 'argv'):
            replay_files(changed)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt = json.loads((bundle/'receipt.json').read_text())
            receipt['archive_sha256'] = '0'*64
            (root/'receipt.json').write_text(json.dumps(receipt))
            (root/'evidence.tar.gz').write_bytes((bundle/'evidence.tar.gz').read_bytes())
            with self.assertRaisesRegex(InvalidEvidence, 'archive changed'):
                replay(root)


if __name__ == '__main__':
    unittest.main()
