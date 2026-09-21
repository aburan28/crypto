import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
import zlib
import evidence


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.content = gzip.compress(b'profile with full instruction costs\n' * 200, mtime=17)
        self.entry = {'path': 'research/test/profile.gz', 'codec': 'gzip-deflate9',
                      'header': self.content[:10].hex(), 'trailer': self.content[-8:].hex(),
                      'bytes': len(self.content), 'mode': 0o644,
                      'sha256': hashlib.sha256(self.content).hexdigest(),
                      'zlib': zlib.ZLIB_RUNTIME_VERSION}

    def bundle(self, entry=None, include_payload=True):
        entry = self.entry if entry is None else entry
        path = self.root / 'test.tar.xz'
        with tarfile.open(path, 'w:xz') as tar:
            index = json.dumps({'format': evidence.FORMAT, 'files': [entry]}).encode()
            info = tarfile.TarInfo('index.json'); info.size = len(index)
            tar.addfile(info, io.BytesIO(index))
            if include_payload:
                payload = gzip.decompress(self.content)
                info = tarfile.TarInfo('payload/' + entry['path']); info.size = len(payload)
                tar.addfile(info, io.BytesIO(payload))
        return path, hashlib.sha256(path.read_bytes()).hexdigest()

    def test_restore_retains_exact_original_gzip_and_is_idempotent(self):
        path, digest = self.bundle()
        destination = self.root / 'restored'
        self.assertEqual(evidence.restore(path, digest, destination)['files'], 1)
        self.assertEqual((destination / self.entry['path']).read_bytes(), self.content)
        self.assertEqual(evidence.restore(path, digest, destination)['status'], 'VERIFIED')

    def test_bad_archive_hash_rejected(self):
        path, _ = self.bundle()
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            evidence.restore(path, '0' * 64, self.root / 'out')

    def test_wrong_reconstruction_rejected(self):
        path, digest = self.bundle({**self.entry, 'sha256': '0' * 64})
        with self.assertRaisesRegex(ValueError, 'reconstruction failed'):
            evidence.restore(path, digest, self.root / 'out')

    def test_different_existing_file_is_preserved(self):
        path, digest = self.bundle(); destination = self.root / 'out'
        existing = destination / self.entry['path']; existing.parent.mkdir(parents=True)
        existing.write_bytes(b'user change')
        with self.assertRaisesRegex(ValueError, 'Refusing to overwrite'):
            evidence.restore(path, digest, destination)
        self.assertEqual(existing.read_bytes(), b'user change')

    def test_path_escape_rejected(self):
        path, digest = self.bundle({**self.entry, 'path': 'research/../../outside'})
        with self.assertRaisesRegex(ValueError, 'Unsafe evidence path'):
            evidence.restore(path, digest, self.root / 'out')

    def test_missing_payload_rejected(self):
        path, digest = self.bundle(include_payload=False)
        with self.assertRaisesRegex(ValueError, 'Incomplete evidence archive'):
            evidence.restore(path, digest, self.root / 'out')

    def multipart(self):
        path, digest = self.bundle()
        data = path.read_bytes()
        middle = len(data) // 2
        row = {'archive': path.name, 'sha256': digest, 'bytes': len(data), 'parts': []}
        for number, content in enumerate((data[:middle], data[middle:]), 1):
            part = self.root / (path.name + '.part' + str(number))
            part.write_bytes(content)
            row['parts'].append({'file': part.name, 'bytes': len(content),
                                 'sha256': hashlib.sha256(content).hexdigest()})
        path.unlink()
        return row

    def test_multipart_retains_exact_original_bytes(self):
        row = self.multipart()
        destination = self.root / 'out'
        result = evidence.restore_bundle(row, self.root, destination)
        self.assertEqual(result, {'archive': row['archive'], 'files': 1,
                                  'status': 'VERIFIED', 'parts': 2})
        self.assertEqual((destination / self.entry['path']).read_bytes(), self.content)

    def test_corrupt_part_rejected_before_extraction(self):
        row = self.multipart()
        part = self.root / row['parts'][0]['file']
        part.write_bytes(b'X' * part.stat().st_size)
        with self.assertRaisesRegex(ValueError, 'part checksum mismatch'):
            evidence.restore_bundle(row, self.root, self.root / 'out')
        self.assertFalse((self.root / 'out').exists())

    def test_missing_part_rejected_before_extraction(self):
        row = self.multipart()
        (self.root / row['parts'][1]['file']).unlink()
        with self.assertRaises(FileNotFoundError):
            evidence.restore_bundle(row, self.root, self.root / 'out')
        self.assertFalse((self.root / 'out').exists())

    def test_reordered_parts_fail_original_archive_hash(self):
        row = self.multipart()
        row['parts'].reverse()
        with self.assertRaisesRegex(ValueError, 'Archive checksum mismatch'):
            evidence.restore_bundle(row, self.root, self.root / 'out')
        self.assertFalse((self.root / 'out').exists())


if __name__ == '__main__':
    unittest.main()
