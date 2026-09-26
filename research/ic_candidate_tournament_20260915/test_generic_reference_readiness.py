"""The five-cell readiness evidence must replay from its frozen archive."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

HARNESS = Path(__file__).parent
RESULT = HARNESS / 'goal_20260924/generic-reference-readiness'


class ReferenceReadinessTests(unittest.TestCase):
    def test_frozen_readiness_archive_replays_all_fifteen_jobs(self):
        manifest = json.loads((RESULT / 'EVIDENCE.json').read_text())
        archive = (RESULT / manifest['archive']).resolve()
        self.assertEqual(hashlib.sha256(archive.read_bytes()).hexdigest(), manifest['archive_sha256'])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with tarfile.open(archive) as source:
                source.extractall(directory, filter='data')
            root = directory / 'ic-generic-reference-readiness'
            for name, expected in manifest['files'].items():
                data = (root / name).read_bytes()
                self.assertEqual(len(data), expected['bytes'])
                self.assertEqual(hashlib.sha256(data).hexdigest(), expected['sha256'])
            output = directory / 'replay.json'
            script = root / 'checkers/goal_20260924/generic-reference-readiness/replay.py'
            process = subprocess.run([sys.executable, str(script), '--evidence', str(root),
                                      '--out', str(output)], capture_output=True, text=True, timeout=120)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            self.assertEqual(output.read_bytes(), (RESULT / 'RESULTS.json').read_bytes())


if __name__ == '__main__':
    unittest.main()
