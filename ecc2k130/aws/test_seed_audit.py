import io
import json
import struct
import unittest

import seed_audit
import seed_registry as sr


class Exists(Exception):
    response = {'Error': {'Code': 'PreconditionFailed'}}


class FakeS3:
    def __init__(self):
        self.data = {}
        self.writes = []

    def put_object(self, Bucket, Key, Body, ContentType, IfNoneMatch):
        assert IfNoneMatch == '*'
        if Key in self.data:
            raise Exists()
        self.data[Key] = Body
        self.writes.append(Key)

    def get_object(self, Bucket, Key):
        return {'Body': io.BytesIO(self.data[Key])}


class BootstrapTests(unittest.TestCase):
    def audit(self):
        return dict(bucket='bucket', records=[dict(version=1, run_id=8000,
                    stream=sr.modal_stream(8000), started=True, checkpoint_floor=10,
                    active_owner='existing-live-owner')])

    def test_readiness_is_written_only_after_every_historical_binding(self):
        s3 = FakeS3()
        seed_audit.initialize(s3, self.audit())
        self.assertEqual(s3.writes, [sr.run_key(8000), sr.CONFIG_KEY])

    def test_rerunning_bootstrap_never_clears_an_owner_or_rewinds_progress(self):
        s3 = FakeS3()
        seed_audit.initialize(s3, self.audit())
        row = json.loads(s3.data[sr.run_key(8000)])
        row.update(active_owner='new-owner', checkpoint_floor=100)
        s3.data[sr.run_key(8000)] = json.dumps(row).encode()
        seed_audit.initialize(s3, self.audit())
        self.assertEqual(json.loads(s3.data[sr.run_key(8000)]), row)

    def test_conflicting_existing_binding_prevents_readiness(self):
        s3 = FakeS3()
        row = dict(run_id=8000, stream='another-stream')
        s3.data[sr.run_key(8000)] = json.dumps(row).encode()
        with self.assertRaisesRegex(ValueError, 'binding differs'):
            seed_audit.initialize(s3, self.audit())
        self.assertNotIn(sr.CONFIG_KEY, s3.data)

    def test_unmapped_legacy_corpus_reads_every_record_not_just_the_first(self):
        records = b''.join(struct.pack('<Q', rid << 48) + b'x' * 24
                           for rid in [4243, 4243, 8000, 9000])
        self.assertEqual(seed_audit.record_run_ids(records), {4243, 8000, 9000})

    def test_v2_corpus_uses_its_header_and_stride(self):
        records = b'ECC2KDP2' + b'\0' * 8 + b''.join(
            struct.pack('<Q', rid << 48) + b'x' * 64 for rid in [12000, 12001])
        self.assertEqual(seed_audit.record_run_ids(records), {12000, 12001})

    def test_truncation_prevents_an_incomplete_seed_inventory(self):
        with self.assertRaisesRegex(ValueError, 'unaligned'):
            seed_audit.record_run_ids(b'x' * 33)


if __name__ == '__main__':
    unittest.main()
