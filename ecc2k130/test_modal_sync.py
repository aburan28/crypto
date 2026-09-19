#!/usr/bin/env python3
import hashlib
import os
import re
import tempfile
import time
import unittest
from pathlib import Path

import modal_sync

ORBIT_KEY_RE = re.compile(
    r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$"
)
CKPT_RE = re.compile(r"^ckpt/(retired/)?slot-(\d+)(?:/([0-9a-f]{64}))?\.ck$")


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

    def test_checkpoint_key_matches_ingest_pattern(self):
        with tempfile.NamedTemporaryFile(delete=False) as fh:
            fh.write(b"checkpoint-bytes")
            path = fh.name
        try:
            key = modal_sync.checkpoint_key(90001, path)
            self.assertTrue(CKPT_RE.match(key), key)
            self.assertEqual(key, "ckpt/slot-90001/%s.ck" % modal_sync.sha256_file(path))
        finally:
            os.unlink(path)

    def test_remote_checkpoint_path(self):
        self.assertEqual(modal_sync.remote_checkpoint(131, 1), "ckpt/curve131-run1.ck")

    def test_run_sh_defaults_a_one_minute_checkpoint(self):
        script = Path(__file__).with_name("run.sh").read_text()
        self.assertIn("CHECKPOINT_EVERY=${CHECKPOINT_EVERY:-60}", script)
        self.assertIn("SYNC_INTERVAL=${SYNC_INTERVAL:-15}", script)
        self.assertIn("--checkpoint-every \"$CHECKPOINT_EVERY\"", script)
        self.assertIn('--watch "$SYNC_INTERVAL"', script)

    def test_checkpoint_header_roundtrip(self):
        blob = modal_sync.pack_checkpoint_header(2, 131, 385024, 16, 1, 1, 6439936)
        self.assertEqual(len(blob), 40)
        parsed = modal_sync.unpack_checkpoint_header(blob)
        self.assertEqual(parsed["threads"], 385024)
        self.assertEqual(parsed["iterBase"], 6439936)
        self.assertEqual(parsed["version"], 2)

    def test_iter_base_from_progress_inverts_the_client_line(self):
        # resume 6439936, then 6160384 walks * 256 steps = 1_577_058_304 pass iters
        self.assertEqual(
            modal_sync.iter_base_from_progress(6439936, 1577058304, 6160384),
            6439936 + 256)

    def test_status_header_file_is_forty_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "curve131-run1.hdr")
            blob = modal_sync.write_checkpoint_header_file(path, dict(
                version=2, m=131, threads=4, batch=16, lanes=1, runId=1, iterBase=32))
            self.assertEqual(os.path.getsize(path), 40)
            with open(path, "rb") as fh:
                self.assertEqual(fh.read(), blob)

    def test_sync_prefers_the_header_sidecar(self):
        with tempfile.TemporaryDirectory() as tmp:
            hdr = os.path.join(tmp, "curve131-run1.hdr")
            modal_sync.write_checkpoint_header_file(hdr, dict(
                version=2, m=131, threads=4, batch=16, lanes=1, runId=1, iterBase=99))
            got = []

            def fake_get(volume, remote, local):
                got.append(remote)
                if remote.endswith(".hdr"):
                    with open(hdr, "rb") as src, open(local, "wb") as dst:
                        dst.write(src.read())
                    return True
                return False

            original = modal_sync.modal_volume_get
            modal_sync.modal_volume_get = fake_get
            try:
                result = modal_sync.sync_checkpoint(
                    None, "bucket", "vol", 131, 1, tmp, dry_run=True)
            finally:
                modal_sync.modal_volume_get = original
            self.assertTrue(result["checkpoint_uploaded"])
            self.assertEqual(got[0], "ckpt/curve131-run1.hdr")
            self.assertTrue(got[0].endswith(".hdr"))

    def test_missing_header_does_not_refetch_the_full_checkpoint_every_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = modal_sync.state_path(131, 1, tmp)
            modal_sync.save_state(state, {
                "checkpoint_sha256": "abc",
                "checkpoint_key": "ckpt/slot-90001/abc.ck",
                "full_checkpoint_check": time.time(),
            })
            got = []

            def fake_get(volume, remote, local):
                got.append(remote)
                return False

            original = modal_sync.modal_volume_get
            modal_sync.modal_volume_get = fake_get
            try:
                result = modal_sync.sync_checkpoint(
                    None, "bucket", "vol", 131, 1, tmp, dry_run=True)
            finally:
                modal_sync.modal_volume_get = original
            self.assertEqual(got, ["ckpt/curve131-run1.hdr"])
            self.assertFalse(result["checkpoint_uploaded"])
            self.assertEqual(result["checkpoint_key"], "ckpt/slot-90001/abc.ck")

    def test_parse_run_ids(self):
        self.assertEqual(modal_sync.parse_run_ids("1,2, 3"), [1, 2, 3])
        self.assertEqual(modal_sync.parse_run_ids(""), [])

    def test_discover_run_ids_from_volume_listing(self):
        class Proc:
            returncode = 0
            stdout = "dp/curve131-run2.bin\ndp/curve131-run1.bin\ndp/curve97-run1.bin\n"
            stderr = ""

        original = modal_sync.subprocess.run
        modal_sync.subprocess.run = lambda *a, **k: Proc()
        try:
            self.assertEqual(modal_sync.discover_run_ids("ecc2k130", 131), [1, 2])
        finally:
            modal_sync.subprocess.run = original


if __name__ == "__main__":
    raise SystemExit(unittest.main())
