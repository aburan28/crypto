"""Count integrity tests for the actual audit helper, without Modal or a GPU."""
import ast
from pathlib import Path
import re
from types import SimpleNamespace
import unittest

from benchreport import benchResult, summarizeSamples


class PackedAuditCountTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1] / 'packed_audit.py'
        helper = next(node for node in ast.parse(path.read_text()).body
                      if isinstance(node, ast.FunctionDef) and node.name == 'checkScalarCounts')
        namespace = dict(re=re, client=SimpleNamespace(benchResult=benchResult))
        exec(compile(ast.Module(body=[helper], type_ignores=[]), str(path), 'exec'), namespace)
        cls.check = staticmethod(namespace['checkScalarCounts'])

    def raw(self, workers=192512, weight=0):
        walks = workers * 32
        expected = walks * 1024 * 32
        return (f'backend cuda-packed131: {workers} threads x 32 slots x 1 lanes = {walks} walks, '
                f'dp weight {weight}, 1024 steps per launch\n'
                f'1.0 s 6000.000 M it/s {expected // 2} iterations 0 dp 0 stored 0 dropped\n'
                f'2.0 s 6000.000 M it/s {expected} iterations 0 dp 0 stored 0 dropped\n'
                'finished: 6000.000 M it/s, 0 distinguished points (0 verified against the reference, 0 dropped)\n')

    def sample(self, raw):
        return benchResult(['fixture'], 0, raw)

    def test_explicit_and_automatic_workers(self):
        for requested in (0, 192512):
            sample = self.sample(self.raw())
            self.assertTrue(self.check(sample, requested, 0))
            self.assertEqual(sample['actualWorkers'], 192512)
            self.assertEqual(sample['expectedIterations'], 201863462912)
            self.assertEqual(sample['reportedIterations'], 201863462912)

    def test_invalid_work_cannot_keep_a_positive_summary(self):
        raw = self.raw()
        first, rest = raw.split('\n', 1)
        variants = {
            'inflated': raw.replace('201863462912 iterations', '403726925824 iterations'),
            'partial': raw.replace('201863462912 iterations', '151397597184 iterations'),
            'ignored workers': self.raw(workers=96256),
            'missing identity': rest,
            'duplicate identity': first + '\n' + raw,
            'wrong batch': raw.replace('32 slots', '16 slots'),
            'wrong steps': raw.replace('1024 steps', '512 steps'),
            'wrong cutoff': self.raw(weight=34),
            'drops': raw.replace('0 dropped', '1 dropped'),
            'nonfinite': raw.replace('finished: 6000.000', 'finished: nan'),
        }
        for label, text in variants.items():
            with self.subTest(label=label):
                sample = self.sample(text)
                self.assertFalse(self.check(sample, 192512, 0))
                self.assertEqual(sample['rate'], 0)
                summary = summarizeSamples([self.sample(raw), sample])
                self.assertFalse(summary['valid'])
                self.assertEqual(summary['rate'], 0)

    def test_collection_uses_its_actual_mode_and_count(self):
        sample = self.sample(self.raw(weight=34))
        self.assertTrue(self.check(sample, 192512, 34))
        self.assertEqual(sample['reportedIterations'], 201863462912)
        self.assertFalse(self.check(self.sample(self.raw()), 192512, 34))


if __name__ == '__main__':
    unittest.main()
