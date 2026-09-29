"""A real retained report exercises the frozen auditor's float-hash failure."""
import json
from pathlib import Path
import tempfile
import unittest

from identity import sha256
from measurement import report_sha256
from oracle import InvalidEvidence
from recover_generic_backend_yield import FROZEN_AUDITOR_SHA256, recover


HERE = Path(__file__).resolve().parent


class FrozenAuditorRepairTests(unittest.TestCase):
    def test_retained_report_requires_the_report_hash_not_the_identity_hash(self):
        source = HERE / 'goal_20260924/generic-public-inputs/final-worker-raw.jsonl'
        with source.open() as stream:
            report = next(json.loads(line)['report'] for line in stream if '"report"' in line)
        self.assertIsInstance(report['elapsed_seconds'], float)
        with self.assertRaises(InvalidEvidence):
            sha256(report)
        self.assertEqual(len(report_sha256(report)), 64)

    def test_changed_archived_auditor_is_rejected_before_import(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            evaluator = bundle / 'tournament/evaluator'
            evaluator.mkdir(parents=True)
            (evaluator / 'generic_backend_yield.py').write_text('unexpected evaluator\n')
            (bundle / 'tournament/contract.json').write_text(json.dumps({
                'qualification_reference_schema': 1, 'seed': 2026092901,
                'evaluator_sha256': {'generic_backend_yield.py': FROZEN_AUDITOR_SHA256}}))
            with self.assertRaisesRegex(ValueError, 'exact frozen'):
                recover(bundle)


if __name__ == '__main__':
    unittest.main()
