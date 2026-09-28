"""Durable mixed-control evidence retains failures and replays without workers."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).parent / 'goal_20260924/generic-adapter-control'
spec = importlib.util.spec_from_file_location('generic_control_archive', HERE / 'audit_archive.py')
archive = importlib.util.module_from_spec(spec)
spec.loader.exec_module(archive)


class GenericControlArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='ic-adapter-archive-test-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = archive.restore(Path(cls.temporary.name))

    def test_archive_hashes_and_export_preserve_all_control_outcomes(self):
        self.assertEqual(archive.control_result(self.root), archive.read(HERE / 'RESULTS.json'))

    @unittest.skipUnless(sys.platform == 'linux', 'strict frozen Linux statistics replay runs in Linux CI')
    def test_frozen_evaluator_reconstructs_all_102_receipts(self):
        self.assertEqual(archive.replay(self.root)['trial_receipts'], 102)


if __name__ == '__main__':
    unittest.main()
