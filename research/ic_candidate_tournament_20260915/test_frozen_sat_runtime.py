"""Historical replay preserves sources without relaxing new execution gates."""
import importlib
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from frozen_sat_runtime import SNAPSHOT, replay, verified_sources
from tournament import read


class FrozenSatRuntimeTests(unittest.TestCase):
    def test_changed_archive_is_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            snapshot = Path(temporary)/'snapshot'
            shutil.copytree(SNAPSHOT, snapshot)
            with (snapshot/'runtime.tar.gz').open('ab') as archive:
                archive.write(b'changed')
            with self.assertRaisesRegex(ValueError, 'archive changed'):
                verified_sources(snapshot)

    def test_unknown_version_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unknown frozen SAT'):
            replay('v3')

    def test_live_execution_still_requires_exact_registered_sources(self):
        for name in ('run_static_sat_full', 'run_static_sat_full_v2'):
            with self.subTest(runner=name):
                runner = importlib.import_module(name)
                with patch.object(runner, 'source_manifest', return_value={}):
                    with self.assertRaisesRegex(ValueError,
                                                'executed source changed'):
                        runner.admit(read(runner.PANEL),
                                     require_local_binary=False)


if __name__ == '__main__':
    unittest.main()
