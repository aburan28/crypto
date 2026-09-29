"""Fresh retained evidence reproduces both registered measurement paths."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location('generic_reference_archive', HERE / 'audit_archive.py')
archive = importlib.util.module_from_spec(spec)
spec.loader.exec_module(archive)


class GenericReferenceArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='ic-generic-reference-archive-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = archive.restore(Path(cls.temporary.name))

    def test_fresh_archive_reproduces_all_exported_runs_and_preserves_scope(self):
        result = archive.exported(self.root)
        self.assertEqual(result, archive.read(HERE / 'RESULTS.json'))
        self.assertEqual(len(result['development_runs']), 990)
        self.assertFalse(result['promotion_eligible'])
        self.assertFalse(result['accepted_reference_binding_changed'])
        self.assertIsNone(result['observer']['legacy_scientific_admission'])

    @unittest.skipUnless(sys.platform == 'linux', 'strict frozen Linux statistics replay runs in Linux CI')
    def test_frozen_evaluators_reconstruct_comparison_and_observer_without_workers(self):
        result = archive.replay(self.root)
        self.assertEqual(result['comparison']['trial_receipts'], 1350)
        self.assertEqual(result['observer']['observer_pairs'], 360)
        self.assertEqual(result['frozen_audits'], 11)
        self.assertEqual(result['distinct_record_ids'], 2550)
        self.assertEqual(result['workers_executed'], 0)


if __name__ == '__main__':
    unittest.main()
