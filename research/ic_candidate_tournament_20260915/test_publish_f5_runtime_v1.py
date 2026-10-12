"""Failure-before-extraction controls for the F4/F5 evidence transport."""
import gzip
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from oracle import InvalidEvidence
from publish_f5_runtime_v1 import publish, replay
from publish_static_sat_v3_control import digest, inventory


class PublishF5RuntimeV1Tests(unittest.TestCase):
    def test_external_invocation_mismatch_cannot_publish_or_extract(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            registration = root / 'registration'
            registration.mkdir()
            (registration / 'execution.json').write_text(json.dumps({'not': 'the registered invocation'}))
            with self.assertRaisesRegex(InvalidEvidence, 'externally frozen'):
                publish(registration, root / 'execution', root / 'audit.json',
                        root / 'published', '0'*64)
            self.assertFalse((root / 'published').exists())

    def test_archive_hash_paths_duplicates_permissions_and_inventory_fail_closed(self):
        cases = [('changed', ['good'], 0o444), ('traversal', ['../escape'], 0o444),
                 ('absolute', ['/escape'], 0o444), ('duplicate', ['same', 'same'], 0o444),
                 ('writable', ['good'], 0o666), ('inventory', ['good'], 0o444)]
        for label, names, mode in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                bundle = root / 'bundle'
                bundle.mkdir()
                buffer = io.BytesIO()
                with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0, filename='') as compressed:
                    with tarfile.open(fileobj=compressed, mode='w|') as tar:
                        for name in names:
                            item = tarfile.TarInfo(name)
                            item.size, item.mode = 4, mode
                            tar.addfile(item, io.BytesIO(b'data'))
                data = buffer.getvalue()
                (bundle / 'evidence.tar.gz').write_bytes(data + (b'changed' if label == 'changed' else b''))
                files = {name: (b'data', mode) for name in names}
                receipt = dict(execution_sha256='a'*64, archive_sha256=digest(data),
                               archive_bytes=len(data), inventory=inventory(files))
                if label == 'inventory':
                    receipt['inventory']['good']['sha256'] = '0'*64
                (bundle / 'receipt.json').write_text(json.dumps(receipt))
                with self.assertRaises(InvalidEvidence):
                    replay(bundle, root / 'output', 'a'*64)
                self.assertFalse((root / 'output').exists())
                self.assertFalse((root / 'escape').exists())


if __name__ == '__main__':
    unittest.main()
