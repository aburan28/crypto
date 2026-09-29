"""The partial result must survive as one verified upload file."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from pack_partial_campaign import pack


class PackPartialCampaignTests(unittest.TestCase):
    def test_preserves_partial_receipt_and_marks_missing_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'campaign'
            receipt = source/'tournament/runs/smoke/one/receipt.json'
            receipt.parent.mkdir(parents=True)
            receipt.write_text('{"status":"TIMEOUT"}\n')
            lock = root/'Cargo.lock'
            lock.write_text('pinned lock\n')
            archive = root/'partial.tar.zst'
            result = pack(source, archive, cargo_lock=lock)
            self.assertFalse(result['capture']['summary_present'])
            self.assertFalse(result['capture']['gate_present'])
            with archive.open('rb') as stream:
                self.assertEqual(result['sha256'], hashlib.file_digest(stream, 'sha256').hexdigest())
            self.assertEqual(result, json.loads((root/'partial.tar.zst.manifest.json').read_text()))
            tar = subprocess.run(['zstd', '-dc', str(archive)], capture_output=True, check=True)
            listing = subprocess.run(['tar', '-tf', '-'], input=tar.stdout,
                                     capture_output=True, check=True, text=False).stdout.decode()
            self.assertIn('campaign/tournament/runs/smoke/one/receipt.json', listing)
            self.assertIn('campaign/workflow-Cargo.lock', listing)
            self.assertIn('campaign/capture.json', listing)


if __name__ == '__main__':
    unittest.main()
