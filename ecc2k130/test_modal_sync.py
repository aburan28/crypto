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
        self.assertIn("VERIFY=${VERIFY:-0}", script)
        self.assertIn("--checkpoint-every \"$CHECKPOINT_EVERY\"", script)
        self.assertIn("--verify $VERIFY", script)
        self.assertIn('--watch "$SYNC_INTERVAL"', script)

    def test_modal_search_defaults_verify_off(self):
        src = Path(__file__).with_name("modal_app.py").read_text()
        self.assertIn("verify: int = 0, checkpoint_every: int = 60, off_campaign: bool = False,\n           cpu_threads: int = -1)", src)
        self.assertIn("packed=False, verify=0,\n              offCampaign=False, cpuThreads=None):", src)

    # --- campaign admission: run ids -------------------------------------

    def test_run_id_range_cannot_reach_an_aws_slot(self):
        # AWS slot s runs as run id s + 1 and slots are 16-bit; the Modal range
        # starts far above anything the fleet has claimed and stays a 5-digit
        # slot under MODAL_SLOT_BASE.
        self.assertGreater(modal_sync.MODAL_RUN_ID_MIN, 4096)
        self.assertEqual(modal_sync.slot_for_run(modal_sync.MODAL_RUN_ID_MAX), 99999)
        with self.assertRaises(ValueError):
            modal_sync.slot_for_run(modal_sync.MODAL_RUN_ID_MAX + 1)

    def test_aws_slot_prefixes_name_every_trace_the_fleet_leaves(self):
        self.assertEqual(modal_sync.aws_slot_prefixes(3),
                         ["ckpt/slot-00002", "ckpt/retired/slot-00002", "dp/slot-00002/"])
        self.assertEqual(modal_sync.aws_slot_prefixes(0), [])

    class BucketWith:
        def __init__(self, keys):
            self.keys = sorted(keys)
            self.asked = []

        def list_objects_v2(self, Bucket, Prefix, MaxKeys=1):
            self.asked.append(Prefix)
            hits = [k for k in self.keys if k.startswith(Prefix)]
            return {"Contents": [{"Key": k} for k in hits[:MaxKeys]]}

    def setUp(self):
        modal_sync._aws_evidence.clear()

    def test_a_run_id_an_aws_slot_has_used_is_refused(self):
        # The 2026-09-20 case: Modal run 3 against AWS slot 2, whose checkpoint
        # is under the legacy pointer name.
        s3 = self.BucketWith(["ckpt/slot-00002.ck", "dp/slot-00002/1789311001-0.bin"])
        reason = modal_sync.admit_run(s3, "bucket", 131, 3)
        self.assertIsNotNone(reason)
        self.assertIn("slot 2", reason)
        self.assertIn("ckpt/slot-00002.ck", reason)

    def test_a_retired_or_orbit_era_slot_is_evidence_too(self):
        s3 = self.BucketWith(["ckpt/retired/slot-00003.ck"])
        self.assertIn("slot 3", modal_sync.admit_run(s3, "bucket", 131, 4))
        modal_sync._aws_evidence.clear()
        s3 = self.BucketWith(["ckpt/slot-00140/%s.ck" % ("ab" * 32)])
        self.assertIn("slot 140", modal_sync.admit_run(s3, "bucket", 131, 141))

    def test_a_neighbouring_slot_is_not_evidence(self):
        # slot-00002 must not match slot-00020: the zero padding is what makes
        # the prefix exact.
        s3 = self.BucketWith(["ckpt/slot-00020.ck", "dp/slot-00021/x-0.bin"])
        reason = modal_sync.admit_run(s3, "bucket", 131, 3)
        self.assertIsNotNone(reason)
        self.assertNotIn("evidence", reason)
        self.assertIn("outside the Modal campaign range", reason)

    def test_a_run_id_in_range_with_no_aws_trace_is_admitted(self):
        s3 = self.BucketWith(["ckpt/slot-00002.ck"])
        self.assertIsNone(modal_sync.admit_run(s3, "bucket", 131, 8000))
        self.assertIsNone(modal_sync.admit_run(s3, "bucket", 131, 9999))

    def test_a_run_id_below_the_range_is_refused_even_with_no_aws_trace_yet(self):
        # An AWS slot could still claim it tomorrow.
        s3 = self.BucketWith([])
        reason = modal_sync.admit_run(s3, "bucket", 131, 4242)
        self.assertIn("outside the Modal campaign range", reason)

    def test_other_curves_are_not_the_campaign(self):
        s3 = self.BucketWith(["ckpt/slot-00002.ck"])
        self.assertIsNone(modal_sync.admit_run(s3, "bucket", 97, 3))

    def test_overrides_are_explicit_and_separate(self):
        s3 = self.BucketWith(["ckpt/slot-00002.ck"])
        lifted = modal_sync.Policy(allow_aws_seed_overlap=True)
        # Overlap lifted, but 3 is still outside the range.
        self.assertIn("outside", modal_sync.admit_run(s3, "bucket", 131, 3, lifted))
        both = modal_sync.Policy(allow_aws_seed_overlap=True, allow_legacy_run_ids=True)
        self.assertIsNone(modal_sync.admit_run(s3, "bucket", 131, 3, both))

    def test_without_a_client_the_range_rule_still_applies(self):
        self.assertIn("outside", modal_sync.admit_run(None, "bucket", 131, 3))
        self.assertIsNone(modal_sync.admit_run(None, "bucket", 131, 8000))

    def test_evidence_is_remembered_and_absence_is_re_asked(self):
        s3 = self.BucketWith(["ckpt/slot-00002.ck"])
        modal_sync.cached_aws_slot_evidence(s3, "bucket", 3, now=1000.0)
        modal_sync.cached_aws_slot_evidence(s3, "bucket", 3, now=5000.0)
        self.assertEqual(len(s3.asked), 1)
        empty = self.BucketWith([])
        modal_sync.cached_aws_slot_evidence(empty, "bucket", 8000, now=1000.0)
        modal_sync.cached_aws_slot_evidence(empty, "bucket", 8000, now=1100.0)
        self.assertEqual(len(empty.asked), 3)
        modal_sync.cached_aws_slot_evidence(empty, "bucket", 8000,
                                            now=1000.0 + modal_sync.AWS_EVIDENCE_TTL_S + 1)
        self.assertEqual(len(empty.asked), 6)

    # --- campaign admission: distinguished-point weight ------------------

    def header(self, iterBase, threads=385024, batch=16, lanes=1):
        return dict(version=2, m=131, threads=threads, batch=batch, lanes=lanes,
                    runId=8000, iterBase=iterBase)

    def test_iterations_multiply_lanes_as_the_client_does(self):
        # walksPerLaunch is threads x BATCH x LANES; a bitsliced header carries
        # LANES > 1 and a packed one carries 1.
        self.assertEqual(modal_sync.iterations_from_header(self.header(10, 2, 4, 1)), 80)
        self.assertEqual(modal_sync.iterations_from_header(self.header(10, 2, 4, 32)), 2560)

    def test_the_theoretical_interval_matches_the_binomial_tail(self):
        self.assertAlmostEqual(modal_sync.theoretical_iter_per_dp_log2(32), 29.01, places=2)
        self.assertAlmostEqual(modal_sync.theoretical_iter_per_dp_log2(34), 25.84, places=2)

    def test_weight_estimates_read_off_the_live_ratios(self):
        # From the bucket on 2026-09-21: the fleet and Modal runs 1-4 at 2^28.41,
        # runs 4243-4245 at 2^25.30 (the documented weight-34 interval), run
        # 4242 at 2^23.58. The listing-inflated readings before re-uploads were
        # deduplicated (2^24.26, 2^21.65) read one weight looser each.
        self.assertEqual(modal_sync.estimate_dp_weight(28.41), 32)
        self.assertEqual(modal_sync.estimate_dp_weight(25.30), 34)
        self.assertEqual(modal_sync.estimate_dp_weight(23.58), 35)
        self.assertEqual(modal_sync.estimate_dp_weight(24.26), 35)
        self.assertIn(modal_sync.estimate_dp_weight(21.65), (36, 37))

    def test_a_campaign_weight_run_passes(self):
        # Run 90002 as measured: 76797696 steps x 6160384 walks over 1324201 records.
        self.assertIsNone(modal_sync.dp_weight_verdict(self.header(76797696), 1324201))

    def test_a_looser_cutoff_is_refused_with_the_weight_named(self):
        # Run 94243's geometry (250000 x 16 walks) at a weight-35 interval.
        verdict = modal_sync.dp_weight_verdict(self.header(8533504, 250000), 1699148)
        self.assertIsNotNone(verdict)
        self.assertIn("looks like weight 35", verdict)
        self.assertIn("--dp-weight 32", verdict)

    def test_too_few_records_is_not_a_verdict(self):
        self.assertIsNone(modal_sync.dp_weight_verdict(self.header(8533504, 250000), 100))
        self.assertIsNone(modal_sync.dp_weight_verdict(None, 10 ** 6))

    def test_sync_once_refuses_before_touching_the_volume(self):
        with tempfile.TemporaryDirectory() as tmp:
            got = []

            def fake_get(volume, remote, local):
                got.append(remote)
                return False

            original = modal_sync.modal_volume_get
            modal_sync.modal_volume_get = fake_get
            try:
                result = modal_sync.sync_once(
                    self.BucketWith(["ckpt/slot-00002.ck"]), "bucket", "vol", 131, 3, tmp,
                    dry_run=True)
            finally:
                modal_sync.modal_volume_get = original
            self.assertIn("slot 2", result["refused"])
            self.assertEqual(result["uploaded_records"], 0)
            self.assertFalse(result["checkpoint_uploaded"])
            self.assertEqual(got, [])

    def test_sync_once_refuses_an_off_weight_run_before_uploading_anything(self):
        with tempfile.TemporaryDirectory() as tmp:
            corpus = os.path.join(tmp, "corpus.bin")
            with open(corpus, "wb") as fh:
                fh.write(b"\0" * 32 * 1699148)
            hdr = os.path.join(tmp, "hdr")
            modal_sync.write_checkpoint_header_file(hdr, self.header(8533504, 250000))

            def fake_get(volume, remote, local):
                src = corpus if remote.endswith(".bin") else hdr if remote.endswith(".hdr") else None
                if src is None:
                    return False
                with open(src, "rb") as s, open(local, "wb") as d:
                    d.write(s.read())
                return True

            uploads = []

            class S3:
                def list_objects_v2(self, **kw):
                    return {"Contents": []}

                def upload_file(self, local, bucket, key):
                    uploads.append(key)

            original = modal_sync.modal_volume_get
            modal_sync.modal_volume_get = fake_get
            try:
                result = modal_sync.sync_once(S3(), "bucket", "vol", 131, 8000,
                                              os.path.join(tmp, "state"))
                self.assertIn("looks like weight 35", result["refused"])
                self.assertEqual(uploads, [])
                lifted = modal_sync.Policy(allow_dp_weight_mismatch=True)
                result = modal_sync.sync_once(S3(), "bucket", "vol", 131, 8000,
                                              os.path.join(tmp, "state"), policy=lifted)
            finally:
                modal_sync.modal_volume_get = original
            self.assertNotIn("refused", result)
            self.assertEqual(result["uploaded_records"], 1699148)
            self.assertTrue(result["checkpoint_uploaded"])
            # Points first, then the header that accounts for them.
            self.assertTrue(uploads[0].startswith("dp/slot-98000/"))
            self.assertTrue(uploads[-1].startswith("ckpt/slot-98000/"))

    def test_main_exposes_the_overrides(self):
        src = Path(__file__).with_name("modal_sync.py").read_text()
        for flag in ("--allow-aws-seed-overlap", "--allow-legacy-run-ids",
                     "--allow-dp-weight-mismatch"):
            self.assertIn(flag, src)

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
