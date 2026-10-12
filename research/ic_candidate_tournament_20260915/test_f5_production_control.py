"""Independent low-width controls for the production-path mathematical auditor."""
import itertools
import copy
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from f5_boolean_control import evaluate
from f5_production_control import decisive_reference, substitute
from f5_production_control import audit
from oracle import InvalidEvidence
from replay_f5_production_control import REGISTRATION, replay, replay_files
from replay_paired_n17_evidence import retained_files


class ProductionControlTests(unittest.TestCase):
    def test_specialization_matches_truth_table_with_cancellation(self):
        rows = [[0, 1, 2, 3, 7], [2, 3], [0, 4, 5], []]
        for variable, value in itertools.product(range(3), (False, True)):
            specialized = substitute(rows, variable, value)
            for assignment in range(8):
                pinned = ((assignment | (1 << variable)) if value
                          else (assignment & ~(1 << variable)))
                self.assertEqual([evaluate(row, pinned) for row in rows],
                                 [evaluate(row, assignment) for row in specialized])

    def test_forces_require_linear_combinations(self):
        _, _, _, contradiction, forced, belongs = decisive_reference([[1, 2], [0, 2]], 3, 3)
        self.assertFalse(contradiction)
        self.assertEqual(sorted(forced), [[0, 1], [0, 2]])
        self.assertTrue(belongs([0, 1]))
        self.assertFalse(belongs([4]))

    def test_nonconstant_consequence_does_not_become_a_force(self):
        _, _, _, contradiction, forced, _ = decisive_reference([[1, 2]], 3, 3)
        self.assertFalse(contradiction)
        self.assertEqual(forced, [])

    def test_constant_refutation_is_detected(self):
        _, _, _, contradiction, _, _ = decisive_reference([[0]], 3, 3)
        self.assertTrue(contradiction)

    def test_nonlinear_decisive_rows_are_sound_on_every_exact_model(self):
        for rows in ([[3, 2], [0, 1]], [[3, 1, 2]], [[7, 1], [0, 2]]):
            models = [a for a in range(8) if all(evaluate(row, a) == 0 for row in rows)]
            self.assertTrue(models)
            _, _, _, contradiction, forced, _ = decisive_reference(rows, 3, 3)
            self.assertFalse(contradiction)
            self.assertTrue(all(evaluate(row, a) == 0 for row in forced for a in models))


class TransportedProductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = Path(__file__).parent/'goal_20260924/f5-production-path-control-20260930/native'
        if not cls.bundle.exists():
            raise unittest.SkipTest('bounded native control not yet archived')
        cls.files = retained_files(cls.bundle)

    def test_actual_native_evidence_replays(self):
        seal = json.loads((self.bundle/'receipt.json').read_text())['archive_sha256']
        result = replay(self.bundle, seal)
        self.assertFalse(result['complete_dlp'])
        self.assertFalse(result['actual_solver_traversal_validated'])
        self.assertEqual(sum(c['unique_systems'] for c in result['controls']), 66)
        self.assertEqual(sum(c['production_matrix_calls'] for c in result['controls']), 132)

    def test_external_archive_seal_is_required(self):
        with self.assertRaises(InvalidEvidence):
            replay(self.bundle, '0'*64)

    def test_changed_outer_archive_bytes_reject(self):
        receipt = json.loads((self.bundle/'receipt.json').read_text())
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            (bundle/'receipt.json').write_text(json.dumps(receipt))
            (bundle/'evidence.tar.gz').write_bytes((self.bundle/'evidence.tar.gz').read_bytes()+b'fault')
            with self.assertRaises(InvalidEvidence):
                replay(bundle, receipt['archive_sha256'])

    def test_binary_kernel_and_raw_export_faults_reject(self):
        for name in ['diagnostic', 'kernel-original.rs', 'native.stdout']:
            with self.subTest(name=name):
                faulty = dict(self.files)
                faulty[name] += b'fault'
                with self.assertRaises(InvalidEvidence):
                    replay_files(faulty)

    def test_changed_source_member_rejects(self):
        faulty = dict(self.files)
        output = io.BytesIO()
        with tarfile.open(fileobj=io.BytesIO(faulty['root-source.tar.gz']), mode='r:gz') as original:
            with tarfile.open(fileobj=output, mode='w:gz') as changed:
                for index, item in enumerate(original):
                    data = original.extractfile(item).read()
                    if index == 0:
                        data += b'changed source member'
                    item.size = len(data)
                    changed.addfile(item, io.BytesIO(data))
        faulty['root-source.tar.gz'] = output.getvalue()
        with self.assertRaises(InvalidEvidence):
            replay_files(faulty)

    def test_changed_native_invocation_rejects(self):
        faulty = dict(self.files)
        process = json.loads(faulty['native-process.json'])
        process['argv'].append('--different-policy')
        faulty['native-process.json'] = json.dumps(process).encode()
        with self.assertRaises(InvalidEvidence):
            replay_files(faulty)

    def test_model_and_specialization_faults_reject_mathematics(self):
        inputs = json.loads(self.files[REGISTRATION+'inputs.json'])
        exported = json.loads(self.files['native.stdout'])
        wrong_model = copy.deepcopy(inputs)
        wrong_model['controls'][0]['models'][0]['assignment'] = '0'
        with self.assertRaises(InvalidEvidence):
            audit(wrong_model, exported)
        wrong_step = copy.deepcopy(exported)
        wrong_step['controls'][0]['paths'][0]['steps'][1][0] = [0]
        with self.assertRaises(InvalidEvidence):
            audit(inputs, wrong_step)

    def test_missing_forced_consequence_rejects_mathematics(self):
        inputs = json.loads(self.files[REGISTRATION+'inputs.json'])
        exported = json.loads(self.files['native.stdout'])
        faulty = copy.deepcopy(exported)
        nodes = faulty['controls'][0]['nodes']
        node = next(n for n in nodes if n['engines'][0]['rows'])
        node['engines'][0]['rows'].pop()
        with self.assertRaises(InvalidEvidence):
            audit(inputs, faulty)


if __name__ == '__main__':
    unittest.main()
