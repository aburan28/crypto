#!/usr/bin/env python3
import hashlib
import os
import re
import tempfile
import unittest

import modal_sync

ORBIT_KEY_RE = re.compile(
    r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$"
)


class ModalSyncTests(unittest.TestCase):
    def test_stream_id_is_stable(self):
        self.assertEqual(len(modal_sync.stream_id(7)), 32)
        self.assertEqual(modal_sync.stream_id(7), modal_sync.stream_id(7))

    def test_slot_mapping(self):
        self.assertEqual(modal_sync.slot_for_run(1), 90001)

    def test_orbit_key_matches_ingest_pattern(self):
        with tempfile.NamedTemporaryFile(delete=False) as fh:
            fh.write(b"x" * 64)
            path = fh.name
        try:
            key = modal_sync.orbit_key(90001, 3, 128, path)
            self.assertTrue(ORBIT_KEY_RE.match(key))
            self.assertEqual(modal_sync.stream_id(3),
                             hashlib.sha256(b"modal-run-3").hexdigest()[:32])
        finally:
            os.unlink(path)

    def test_state_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = modal_sync.state_path(131, 1, tmp)
            modal_sync.save_state(path, {"offset": 320, "uploaded_records": 10})
            self.assertEqual(modal_sync.load_state(path)["offset"], 320)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
