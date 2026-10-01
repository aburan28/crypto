"""Failure-before-extraction controls for the F4/F5 evidence transport."""
import gzip
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from f5_runtime_receipts_v2 import publication_process
from identity import sha256, write_immutable

from oracle import InvalidEvidence
from publish_f5_runtime_v2 import publish, replay
from publish_static_sat_v3_control import digest, inventory


class PublishF5RuntimeV2Tests(unittest.TestCase):
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

    def test_raw_diagnostics_publish_transport_without_float_writer_failure(self):
        # Isolate publication from scientific admission. The actual source/math
        # auditor is covered separately; this control makes no solver claim.
        retained = Path(__file__).resolve().parent/'goal_20260924/f5-source-bound-runtime-v1/results-20260930/NATIVE-PARSE-REPLAY.json'
        raw = json.loads(retained.read_text())['native']
        result = dict(schema_version=2, status='PUBLICATION_CONTROL_ONLY',
            native_pipeline=publication_process(raw), complete_ic_admitted=False,
            online_speedup=None, promotion_eligible=False)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            registration, execution = root/'registration', root/'execution'
            registration.mkdir()
            execution.mkdir()
            spec = {'publication_control_only': True}
            write_immutable(registration/'execution.json', spec)
            raw_bytes = json.dumps(raw, sort_keys=True).encode()
            (execution/'pipeline.metrics.json').write_bytes(raw_bytes)
            audit_file = root/'audit.json'
            write_immutable(audit_file, result)
            with patch('publish_f5_runtime_v2.audit', return_value=result):
                self.assertEqual(publish(registration, execution, audit_file, root/'bundle', sha256(spec)), result)
                self.assertEqual(replay(root/'bundle', root/'transport', sha256(spec)), result)
            self.assertEqual((root/'transport/execution/pipeline.metrics.json').read_bytes(), raw_bytes)
            self.assertEqual(json.loads((root/'bundle/AUDIT.json').read_text()), result)

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
