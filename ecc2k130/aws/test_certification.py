"""Local certification gates: real collision resolution plus injected failures.

No AWS credentials, services, GPU, or network are used by this suite.
"""
import copy
import fcntl
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import merge
import protocol
import worker

ROOT = Path(__file__).resolve().parents[1]
CLIENT = ROOT / "ecc2k130-cpu"
FIXTURE = ROOT / "build/test-production"
RECORD = struct.Struct("<4Q")


class Certification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(subprocess.check_output([str(FIXTURE), "--fixture"], text=True))
        cls.binaryHash = protocol.sha256File(CLIENT)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ecc-cert-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "corpus"
        self.source.mkdir()
        self.work = self.root / "merge"
        self.work.mkdir()
        self.config = dict(storageProtocol=protocol.PROTOCOL, curve=41,
                           dpWeight=self.fixture["weight"], maxIters=100000,
                           packed=False, workers=1, batch=32, extraArgs=[],
                           binarySha256=self.binaryHash, hostBinarySha256=self.binaryHash, sourceSha256="a" * 64,
                           steps=4, checkpointEvery=1, uploadEvery=1, verify=0,
                           restartHours=0, loadMax=100000)
        self.campaign = protocol.campaignContract(self.config)
        self.configPath = self.root / "campaign.json"
        protocol.atomicJson(self.configPath, self.config)

    def chunk(self, name, seeds, committed=True):
        path = self.source / name
        path.write_bytes(b"".join(RECORD.pack(seed, *self.fixture["key"]) for seed in seeds))
        if committed:
            protocol.atomicJson(str(path) + ".json", protocol.envelope(path, self.campaign, "dp"))
        return path

    def callMerge(self, *extra):
        return subprocess.run([sys.executable, str(ROOT / "aws/merge.py"), "--work", str(self.work),
                               "--local", str(self.source), "--campaign", str(self.configPath),
                               "--client", str(CLIENT), "--buckets", "16", *extra],
                              text=True, capture_output=True, timeout=30)

    def client(self, *extra):
        return subprocess.run([str(CLIENT), "--curve", "83", "--threads", "1",
                               "--steps", "1", "--launches", "1", "--verify", "0", *extra],
                              text=True, capture_output=True, timeout=30)

    def fakeWorker(self):
        w = worker.Worker.__new__(worker.Worker)
        w.work = str(self.root / "worker")
        os.mkdir(w.work)
        w.store = worker.LocalStore(str(self.root / "store"))
        os.mkdir(w.store.root)
        w.slots = worker.LocalSlots(str(self.root / "slots.json"))
        w.owner = "test:gpu0:unique-epoch"
        w.statePath = os.path.join(w.work, "state.json")
        w.state = {"dpOffset": 0, "ckptIter": -1}
        w.contract = self.campaign
        w.leaseLost = False
        w.stopping = False
        w.streamId = "test-stream"
        w.lastBeatSuccess = time.monotonic()
        slot = w.slots.claim(w.owner, {"campaignId": self.campaign["id"]})
        Path(w.dpPath).write_bytes(RECORD.pack(self.fixture["seedA"], *self.fixture["key"]))
        Path(w.ckptPath).write_bytes(struct.pack("<8s6IQ", b"ECC2K130", 1, 41, 1, 32, 64, 1, 99))
        return w, slot

    def test_exact_known_cross_corpus_collision_and_idempotence(self):
        a, b = self.fixture["seedA"], self.fixture["seedB"]
        self.chunk("worker-a.bin", [a, a])
        self.chunk("worker-b.bin", [b])
        # Each worker alone holds only one seed, so this cannot pass via a
        # worker-local solve or a newly launched random walk.
        p = self.callMerge("--detect-only")
        self.assertEqual(p.returncode, 0, p.stderr)
        s = json.loads(p.stdout)
        self.assertEqual((s["added"], s["corpus"], s["collisions"]), (3, 2, 1))
        p = self.callMerge()
        self.assertEqual(p.returncode, 0, p.stderr)
        s = json.loads(p.stdout)
        self.assertEqual(s["added"], 0)
        self.assertEqual(s["solution"]["k"], self.fixture["k"])
        self.assertTrue(s["solution"]["verified"])
        p = self.callMerge()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(json.loads(p.stdout)["k"], self.fixture["k"])

    def test_chunk_tampering_rejected_before_offset_commit(self):
        p = self.chunk("a.bin", [self.fixture["seedA"]])
        data = bytearray(p.read_bytes()); data[0] ^= 1; p.write_bytes(data)
        self.assertNotEqual(self.callMerge().returncode, 0)
        self.assertFalse((self.work / "state.json").exists())

    def test_foreign_campaign_manifest_rejected(self):
        p = self.chunk("a.bin", [self.fixture["seedA"]])
        meta = json.loads(Path(str(p) + ".json").read_text())
        meta["campaignId"] = "b" * 64
        protocol.atomicJson(str(p) + ".json", meta)
        self.assertNotEqual(self.callMerge().returncode, 0)

    def test_uncommitted_chunk_not_ingested_then_committed(self):
        p = self.chunk("a.bin", [self.fixture["seedA"]], committed=False)
        r = self.callMerge("--detect-only")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(json.loads(r.stdout)["added"], 0)
        protocol.atomicJson(str(p) + ".json", protocol.envelope(p, self.campaign, "dp"))
        r = self.callMerge("--detect-only")
        self.assertEqual(json.loads(r.stdout)["added"], 1)

    def test_campaign_change_rejected_even_after_solved(self):
        self.chunk("a.bin", [self.fixture["seedA"], self.fixture["seedB"]])
        self.assertEqual(self.callMerge().returncode, 0)
        self.config["dpWeight"] += 1
        protocol.atomicJson(self.configPath, self.config)
        self.assertNotEqual(self.callMerge().returncode, 0)

    def test_merge_exclusive_lock(self):
        with open(self.work / "merge.lock", "a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.assertNotEqual(self.callMerge().returncode, 0)

    def test_bad_bucket_counts_rejected(self):
        for count in (0, -1, 3):
            with self.subTest(count=count):
                self.assertNotEqual(self.callMerge("--buckets", str(count)).returncode, 0)

    def test_crash_after_bucket_append_before_state_is_idempotent(self):
        self.chunk("a.bin", [self.fixture["seedA"], self.fixture["seedB"]])
        state = merge.loadState(str(self.work / "state.json")); state["buckets"] = 16
        before = copy.deepcopy(state)
        merge.ingest(state, str(self.source), str(self.work), self.campaign)
        # Crash: discard updated offsets, retain durable bucket bytes.
        merge.ingest(before, str(self.source), str(self.work), self.campaign)
        found, total = merge.detect(before, str(self.work))
        self.assertEqual((len(found), total), (1, 2))

    def test_torn_bucket_fails_closed(self):
        self.chunk("a.bin", [self.fixture["seedA"]])
        state = merge.loadState(str(self.work / "state.json")); state["buckets"] = 16
        merge.ingest(state, str(self.source), str(self.work), self.campaign)
        path = next((self.work / "buckets").glob("*.bin"))
        with path.open("ab") as fh:
            fh.write(b"torn")
        with self.assertRaisesRegex(ValueError, "torn bucket"):
            merge.detect(state, str(self.work))

    def test_committed_bucket_bit_corruption_detected(self):
        self.chunk("a.bin", [self.fixture["seedA"]])
        self.assertEqual(self.callMerge("--detect-only").returncode, 0)
        path = next((self.work / "buckets").glob("*.bin"))
        data = bytearray(path.read_bytes()); data[10] ^= 1; path.write_bytes(data)
        p = self.callMerge("--detect-only")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("bucket integrity failure", p.stderr)

    def test_unpinned_host_solver_rejected(self):
        config = dict(self.config, hostBinarySha256="b" * 64)
        protocol.atomicJson(self.configPath, config)
        p = self.callMerge()
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("solver binary hash", p.stderr)

    def test_input_truncation_not_silently_skipped(self):
        p = self.chunk("a.bin", [self.fixture["seedA"]])
        state = merge.loadState(str(self.work / "state.json"))
        merge.ingest(state, str(self.source), str(self.work))
        p.write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "shrank"):
            merge.ingest(state, str(self.source), str(self.work))

    def test_timed_out_solve_remains_retryable(self):
        self.chunk("a.bin", [self.fixture["seedA"], self.fixture["seedB"]])
        self.assertEqual(self.callMerge("--detect-only").returncode, 0)
        p = self.callMerge("--solve-timeout", "0.000001")
        self.assertEqual(p.returncode, 2, p.stderr)
        p = self.callMerge()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(json.loads(p.stdout)["solution"]["k"], self.fixture["k"])

    def test_write_open_failure_is_fatal(self):
        self.assertEqual(self.client("--dp-file", str(self.root / "missing/dp.bin")).returncode, 8)

    def test_device_file_not_accepted_as_corpus(self):
        self.assertEqual(self.client("--dp-file", "/dev/null").returncode, 8)

    def test_partial_raw_record_rejected_without_modification(self):
        p = self.root / "partial.bin"; p.write_bytes(b"broken")
        self.assertEqual(self.client("--dp-file", str(p)).returncode, 8)
        self.assertEqual(p.read_bytes(), b"broken")

    def test_report_overflow_never_checkpoints(self):
        ck = self.root / "state.ck"
        p = self.client("--dp-weight", "83", "--dp-cap", "1", "--checkpoint", str(ck))
        self.assertEqual(p.returncode, 7, p.stdout + p.stderr)
        self.assertFalse(ck.exists())

    def test_corpus_single_writer_lock(self):
        p = self.root / "dp.bin"
        with p.open("ab") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.assertEqual(self.client("--dp-file", str(p)).returncode, 8)

    def test_invalid_parameters_rejected(self):
        for args in (("--run-id", "65536"), ("--steps", "0"), ("--dp-cap", "0"),
                     ("--dp-weight", "132"), ("--checkpoint-every", "nan")):
            with self.subTest(args=args):
                self.assertEqual(self.client(*args).returncode, 1)

    def test_contract_rejects_overrides_and_missing_hash(self):
        for key, value in (("extraArgs", ["--dp-weight", "3"]), ("binarySha256", "")):
            config = dict(self.config, **{key: value})
            with self.assertRaises(ValueError):
                protocol.campaignContract(config)

    def test_contract_binds_every_semantic_field(self):
        for key, value in (("dpWeight", 14), ("curve", 83), ("maxIters", 99),
                           ("binarySha256", "b" * 64), ("sourceSha256", "b" * 64),
                           ("workers", 2), ("batch", 16), ("walk", "table")):
            with self.subTest(field=key):
                self.assertNotEqual(protocol.campaignContract(dict(self.config, **{key: value}))["id"], self.campaign["id"])
        # The default is the walk every existing corpus was collected with.
        self.assertEqual(protocol.campaignContract(dict(self.config, walk="sigma"))["id"], self.campaign["id"])
        with self.assertRaises(ValueError):
            protocol.campaignContract(dict(self.config, walk="random"))

    def test_legacy_directory_cannot_be_silently_certified(self):
        (self.work / "dp.bin").write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "legacy"):
            protocol.bindDirectory(self.work, self.campaign)

    def test_snapshot_publication_order_and_integrity(self):
        w, slot = self.fakeWorker()
        published = []
        put = w.store.put
        def capture(src, key):
            published.append(key); return put(src, key)
        with patch.object(w.store, "put", side_effect=capture):
            w.uploadCycle(slot)
        self.assertTrue(published[0].endswith(".bin"))
        self.assertTrue(published[1].endswith(".bin.json"))
        self.assertTrue(published[2].endswith(".ck"))
        self.assertTrue(published[3].endswith(".ck.json"))
        record = w.slots.scan()[0]
        self.assertEqual(record["checkpointKey"], published[2])
        self.assertEqual(record["ckptIter"], 99)
        saved = Path(w.store.root) / published[2]
        protocol.verifyEnvelope(saved, json.loads(Path(str(saved) + ".json").read_text()), self.campaign, "checkpoint")

    def test_failed_dp_upload_cannot_publish_checkpoint(self):
        w, slot = self.fakeWorker()
        with patch.object(w.store, "put", side_effect=OSError("injected upload failure")):
            with self.assertRaises(OSError):
                w.uploadCycle(slot)
        self.assertNotIn("checkpointKey", w.slots.scan()[0])
        self.assertEqual(w.state["dpOffset"], 0)

    def test_failed_manifest_upload_does_not_commit_offset(self):
        w, slot = self.fakeWorker()
        put = w.store.put
        def failManifest(src, key):
            if key.endswith(".json"):
                raise OSError("manifest upload failure")
            return put(src, key)
        with patch.object(w.store, "put", side_effect=failManifest):
            with self.assertRaises(OSError):
                w.uploadCycle(slot)
        self.assertEqual(w.state["dpOffset"], 0)
        self.assertNotIn("checkpointKey", w.slots.scan()[0])
        w.uploadCycle(slot)
        self.assertEqual(w.state["dpOffset"], 32)

    def test_lease_loss_before_upload_prevents_publication(self):
        w, slot = self.fakeWorker()
        w.slots.release(slot, w.owner)
        with patch.object(w.store, "put") as put:
            with self.assertRaisesRegex(RuntimeError, "lease lost"):
                w.uploadCycle(slot)
            put.assert_not_called()

    def test_lease_loss_during_checkpoint_upload_cannot_commit_pointer(self):
        w, slot = self.fakeWorker()
        put = w.store.put
        def loseLease(src, key):
            result = put(src, key)
            if key.endswith(".ck.json"):
                w.slots.release(slot, w.owner)
                w.slots.claim("replacement", {"campaignId": self.campaign["id"]})
            return result
        with patch.object(w.store, "put", side_effect=loseLease):
            with self.assertRaisesRegex(RuntimeError, "pointer not committed"):
                w.uploadCycle(slot)
        record = w.slots.scan()[0]
        self.assertEqual(record["owner"], "replacement")
        self.assertNotIn("checkpointKey", record)

    def test_rotation_crash_can_only_duplicate_not_skip(self):
        w, slot = self.fakeWorker()
        w.uploadCycle(slot)
        unlink = worker.os.remove
        def crash(path):
            if path == w.dpPath:
                raise OSError("injected crash before unlink")
            return unlink(path)
        with patch.object(worker.os, "remove", side_effect=crash):
            with self.assertRaises(OSError):
                w.rotateDpFile()
        self.assertEqual(json.loads(Path(w.statePath).read_text())["dpOffset"], 0)
        self.assertEqual(Path(w.dpPath).stat().st_size, 32)

    def test_checkpoint_manifest_corruption_rejected(self):
        w, slot = self.fakeWorker()
        w.uploadCycle(slot)
        saved = Path(w.store.root) / w.slots.scan()[0]["checkpointKey"]
        meta = json.loads(Path(str(saved) + ".json").read_text())
        data = bytearray(saved.read_bytes()); data[-1] ^= 1; saved.write_bytes(data)
        with self.assertRaisesRegex(ValueError, "manifest mismatch"):
            protocol.verifyEnvelope(saved, meta, self.campaign, "checkpoint")

    def test_expired_lease_cannot_be_revived_by_old_owner(self):
        w, slot = self.fakeWorker()
        with patch.object(worker.time, "time", return_value=time.time() + worker.LEASE_SECONDS + 1):
            self.assertFalse(w.slots.heartbeat(slot, w.owner, {}))

    def test_strict_worker_resume_on_replacement_host(self):
        store = self.root / "resume-store"; store.mkdir()
        config = dict(self.config, curve=83, dpWeight=24, steps=8)
        protocol.atomicJson(store / "campaign.json", config)
        program = ("import worker; worker.instanceId=lambda:'local-cert'; "
                   "worker.gpuName=lambda g:'cpu'; raise SystemExit(worker.Worker().run())")
        lastIteration = -1
        for attempt in range(2):
            root = self.root / ("replacement-%d" % attempt); root.mkdir()
            env = dict(os.environ, ECC_ROOT=str(root), ECC_LOCAL_STORE=str(store),
                       ECC_CLIENT=str(CLIENT), ECC_GPU="0")
            env.pop("ECC_ALLOW_LEGACY_STORAGE", None)
            log = self.root / ("resume-%d.log" % attempt)
            with log.open("w") as out:
                p = subprocess.Popen([sys.executable, "-c", program], cwd=ROOT / "aws", env=env,
                                     text=True, stdout=out, stderr=subprocess.STDOUT)
                deadline = time.monotonic() + 25
                try:
                    while time.monotonic() < deadline and p.poll() is None:
                        statePath = store / "slots.json"
                        rows = json.loads(statePath.read_text()) if statePath.exists() else {}
                        row = rows.get("0", {})
                        if row.get("ckptIter", -1) > lastIteration and row.get("checkpointKey"):
                            break
                        time.sleep(0.05)
                    else:
                        self.fail("replacement did not publish checkpoint: " + log.read_text())
                    p.terminate(); p.wait(timeout=15)
                finally:
                    if p.poll() is None:
                        p.kill(); p.wait()
            self.assertEqual(p.returncode, 0, log.read_text())
            row = json.loads((store / "slots.json").read_text())["0"]
            self.assertGreater(row["ckptIter"], lastIteration)
            lastIteration = row["ckptIter"]
            if attempt:
                self.assertIn("downloaded checkpoint", log.read_text())
                self.assertIn("resumed from", log.read_text())

    def test_s3_access_denied_is_not_missing_checkpoint(self):
        response = subprocess.CompletedProcess([], 1, "", "AccessDenied")
        with patch.object(worker.subprocess, "run", return_value=response):
            with self.assertRaises(RuntimeError):
                worker.S3Store("unused-test-bucket").exists("checkpoint")

    def test_s3_claim_race_does_not_mutate_scanned_records(self):
        slots = worker.S3Slots("unused-test-bucket")
        rows = [{"slot": 0, "_etag": "etag", "state": "idle", "leaseUntil": 0}]
        with patch.object(slots, "scan", return_value=rows), patch.object(slots, "_put", side_effect=[False, True]):
            self.assertEqual(slots.claim("new-owner", {}), 1)
        self.assertEqual(rows[0]["slot"], 0)

    def test_worker_strict_known_log_end_to_end(self):
        # Real Worker + real client + local stand-in for S3/slot registry.
        root, store = self.root / "real-worker", self.root / "real-store"
        root.mkdir(); store.mkdir()
        protocol.atomicJson(store / "campaign.json", self.config)
        env = dict(os.environ, ECC_ROOT=str(root), ECC_LOCAL_STORE=str(store),
                   ECC_CLIENT=str(CLIENT), ECC_GPU="0")
        env.pop("ECC_ALLOW_LEGACY_STORAGE", None)
        # Avoid EC2 metadata and nvidia-smi lookups: the test is strictly offline.
        program = ("import worker; worker.instanceId=lambda:'local-cert'; "
                   "worker.gpuName=lambda g:'cpu'; raise SystemExit(worker.Worker().run())")
        p = subprocess.run([sys.executable, "-c", program], cwd=ROOT / "aws", env=env,
                           text=True, capture_output=True, timeout=45)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        solved = json.loads((store / "solution.json").read_text())
        self.assertIn(self.fixture["k"], solved["line"])
        for path in (store / "dp").rglob("*.bin"):
            protocol.verifyEnvelope(path, json.loads(Path(str(path) + ".json").read_text()), self.campaign, "dp")
        rows = json.loads((store / "slots.json").read_text())
        self.assertEqual(rows["0"]["state"], "solved")


if __name__ == "__main__":
    unittest.main()
