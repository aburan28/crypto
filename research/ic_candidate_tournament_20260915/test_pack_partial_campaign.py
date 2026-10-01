"""The archive survives a concurrent-read warning with a failed claim gate."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
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
            self.assertEqual(result['pack_status'], 'ARCHIVE_READ_COMPLETE')
            self.assertEqual(result['tar_exit_code'], 0)
            self.assertEqual(result['tar_stderr_bytes'], 0)
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

    def test_changed_directory_retains_censored_archive_and_warning(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'campaign'
            source.mkdir()
            archive = root/'partial.tar.zst'
            stream = io.BytesIO()
            with tarfile.open(fileobj=stream, mode='w') as tar:
                body = b'{"status":"TIMEOUT"}\n'
                member = tarfile.TarInfo('campaign/receipt.json')
                member.size = len(body)
                tar.addfile(member, io.BytesIO(body))
            payload = stream.getvalue()

            class ChangedTar:
                def __init__(self, output):
                    self.stdout = output

                def wait(self):
                    return 1

            def changed_directory(*_args, **kwargs):
                kwargs['stderr'].write(b'tar: campaign: file changed as we read it\n')
                output = tempfile.TemporaryFile()
                output.write(payload)
                output.seek(0)
                return ChangedTar(output)

            result = pack(source, archive, tar_factory=changed_directory)
            self.assertEqual(result['pack_status'], 'ARCHIVE_READ_UNVERIFIED')
            self.assertEqual(result['tar_exit_code'], 1)
            self.assertTrue(archive.is_file())
            self.assertGreater(result['tar_stderr_bytes'], 0)
            self.assertTrue((root/result['tar_stderr_file']).is_file())
            self.assertEqual(result, json.loads(
                (root/'partial.tar.zst.manifest.json').read_text()))
            restored = subprocess.run(['zstd', '-dc', str(archive)],
                                      capture_output=True, check=True).stdout
            with tarfile.open(fileobj=io.BytesIO(restored), mode='r:') as tar:
                self.assertEqual(tar.extractfile('campaign/receipt.json').read(),
                                 b'{"status":"TIMEOUT"}\n')


if __name__ == '__main__':
    unittest.main()
