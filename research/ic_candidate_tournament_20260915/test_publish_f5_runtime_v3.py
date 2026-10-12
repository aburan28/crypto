"""Publication boundary controls; no native worker or mathematical admission."""
import gzip
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from oracle import InvalidEvidence
from publish_f5_runtime_v3 import publish, replay
from publish_static_sat_v3_control import digest, inventory


class PublishF5V3Tests(unittest.TestCase):
    def test_wrong_invocation_cannot_read_or_publish_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);reg=root/'registration';reg.mkdir()
            (reg/'execution.json').write_text('{}')
            with self.assertRaisesRegex(InvalidEvidence,'externally frozen'):
                publish(reg,root/'absent-execution',root/'absent-audit',
                        root/'absent-prior',root/'output','0'*64)
            self.assertFalse((root/'output').exists())

    def test_changed_hash_unsafe_members_modes_and_inventory_cannot_extract(self):
        cases=[('hash',['safe'],0o444),('traversal',['../escape'],0o444),
               ('absolute',['/escape'],0o444),('duplicate',['same','same'],0o444),
               ('writable',['safe'],0o666),('inventory',['safe'],0o444)]
        for label,names,mode in cases:
            with self.subTest(label=label),tempfile.TemporaryDirectory() as temporary:
                root=Path(temporary);bundle=root/'bundle';bundle.mkdir()
                stream=io.BytesIO()
                with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0,filename='') as compressed:
                    with tarfile.open(fileobj=compressed,mode='w|') as tar:
                        for name in names:
                            member=tarfile.TarInfo(name);member.size=4;member.mode=mode
                            tar.addfile(member,io.BytesIO(b'data'))
                data=stream.getvalue()
                (bundle/'evidence.tar.gz').write_bytes(data+(b'changed' if label=='hash' else b''))
                inv=inventory({name:(b'data',mode) for name in names})
                if label=='inventory': inv['safe']['sha256']='0'*64
                receipt={'schema_version':3,'execution_sha256':'a'*64,
                         'archive_sha256':digest(data),'archive_bytes':len(data),'inventory':inv}
                (bundle/'receipt.json').write_text(json.dumps(receipt))
                with self.assertRaises(InvalidEvidence): replay(bundle,root/'output','a'*64)
                self.assertFalse((root/'output').exists())
                self.assertFalse((root/'escape').exists())


if __name__=='__main__':
    unittest.main()
