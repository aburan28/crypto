"""Points wait on disk when the store will not take them.

The five-hour outage of 2026-09-17 was the ingest refusing to read the
objects, not the workers failing to write them, but it put the question:
until now an upload that kept failing left the records in dp.bin and nowhere
else, and claimSlot and rotateDpFile are both entitled to delete that.  The
spool makes surviving a failed upload deliberate rather than lucky, for
process and transport failures only -- it is local disk, so a Spot
reclamation still takes whatever is in it.

No network, no AWS: a LocalStore in a temporary directory stands in for S3,
and a store that raises on put stands in for S3 being down.

No type hints, camelCase identifiers (project convention).
"""
import contextlib
import io
import os
import shutil
import struct
import tempfile
import unittest
from unittest import mock

import worker
from protocol import PROTOCOL, campaignContract, sha256File, verifyEnvelope
from worker import (CKPT_MAGIC, DP_HEADER_BYTES, DP_MAGIC_V2, RECORD_BYTES,
                    RECORD_BYTES_V2, Worker, readJson)


def record(i):
    """One 32-byte distinguished point record, distinguishable from the rest."""
    return ("%032d" % i).encode()


def recordV2(i):
    """One 72-byte v2 record: the same point, plus the witness."""
    return ("%072d" % i).encode()


def v2Header():
    return DP_MAGIC_V2 + struct.pack("<II", 2, RECORD_BYTES_V2)


def storedRecordsV2(root):
    """Every record the store holds, framing each object by its own magic.

    The merge frames each delta standalone, so this reads them the same way.
    An object that does not announce v2 is read as v1, which is exactly the
    mis-framing a missing per-delta header would cause.
    """
    out = []
    for dirpath, _dirs, names in os.walk(os.path.join(root, "dp")):
        for name in sorted(names):
            if not name.endswith(".bin"):
                continue
            with open(os.path.join(dirpath, name), "rb") as fh:
                data = fh.read()
            if data[:len(DP_MAGIC_V2)] == DP_MAGIC_V2:
                body, stride = data[DP_HEADER_BYTES:], RECORD_BYTES_V2
            else:
                body, stride = data, RECORD_BYTES
            out += [body[i:i + stride] for i in range(0, len(body), stride)]
    return out


def quiet(fn, *args):
    """Run fn, returning (log output, whatever it raised or None)."""
    out = io.StringIO()
    err = None
    with contextlib.redirect_stdout(out):
        try:
            fn(*args)
        except Exception as e:
            err = e
    return out.getvalue(), err


def storedRecords(root):
    """Every dp record the store holds, in no particular order."""
    out = []
    dp = os.path.join(root, "dp")
    for dirpath, _dirs, names in os.walk(dp):
        for name in sorted(names):
            if not name.endswith(".bin"):
                continue
            with open(os.path.join(dirpath, name), "rb") as fh:
                data = fh.read()
            out += [data[i:i + RECORD_BYTES] for i in range(0, len(data), RECORD_BYTES)]
    return out


def storedKeys(root):
    out = []
    for dirpath, _dirs, names in os.walk(os.path.join(root, "dp")):
        for name in sorted(names):
            out.append(os.path.relpath(os.path.join(dirpath, name), root))
    return sorted(out)


class BrokenStore:
    """S3 is down: lookups still answer, puts do not."""

    def __init__(self, inner):
        self.inner = inner
        self.puts = 0

    def exists(self, key):
        return self.inner.exists(key)

    def get(self, key, dest):
        return self.inner.get(key, dest)

    def put(self, src, key):
        self.puts += 1
        raise RuntimeError("store unreachable")


class CountingStore:
    """A healthy store that remembers which keys were written through it."""

    def __init__(self, inner):
        self.inner = inner
        self.puts = []

    def exists(self, key):
        return self.inner.exists(key)

    def get(self, key, dest):
        return self.inner.get(key, dest)

    def put(self, src, key):
        self.puts.append(key)
        self.inner.put(src, key)


class SpoolCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = os.path.join(self.tmp.name, "root")
        self.storeRoot = os.path.join(self.tmp.name, "store")
        os.makedirs(self.root)
        # ECC_DEVICE with the name, not the name alone: a slot record whose
        # gpuName is cpu-shaped is refused to a claimant that did not declare
        # itself a CPU worker (idleSlotClaimable), so a name-only environment
        # is a GPU worker looking at a CPU slot -- a pair bootstrap-cpu.sh
        # never produces, and one that made this fixture skip its own slot.
        env = mock.patch.dict(os.environ, {"ECC_ROOT": self.root,
                                           "ECC_LOCAL_STORE": self.storeRoot,
                                           "ECC_GPU": "0",
                                           "ECC_DEVICE": "cpu",
                                           "ECC_THREADS": "1",
                                           "ECC_DEVICE_NAME": "cpu",
                                           "ECC_INSTANCE_TYPE": "local"})
        env.start()
        self.addCleanup(env.stop)
        # instanceId() would ask the EC2 metadata service.
        ident = mock.patch.object(worker, "instanceId", lambda: "test-instance")
        ident.start()
        self.addCleanup(ident.stop)
        self.written = 0

    def makeWorker(self):
        w = Worker()
        self.addCleanup(w.workLock.close)
        w.cfg = {}
        return w

    def claimed(self):
        """A worker holding a slot, with a healthy store."""
        w = self.makeWorker()
        _, err = quiet(w.claimSlot)
        self.assertIsNone(err)
        return w, int(w.state["slot"])

    def appendDp(self, w, count):
        with open(w.dpPath, "ab") as fh:
            for i in range(self.written, self.written + count):
                fh.write(record(i))
        self.written += count

    def failedCycle(self, w, slot):
        """One upload cycle against a store that refuses every put."""
        healthy = w.store
        w.store = BrokenStore(healthy)
        out, err = quiet(w.uploadCycle, slot)
        w.store = healthy
        self.assertIsNotNone(err, "a store that raises on put must fail the cycle")
        self.lastFailure = err
        return out


class FailedUpload(SpoolCase):
    def test_a_second_failed_cycle_does_not_stack_a_superset(self):
        """A cut starts at dpOffset, which a failed upload leaves alone.

        Cutting again would spool [offset, more) beside the [offset, less)
        that just failed -- the same records twice, once per cycle, for as
        long as the outage lasts.  The later records are not lost by
        refusing: they are still in dp.bin, which is not rotated while the
        spool holds anything.
        """
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        first = [e["name"] for e in w.spoolEntries()]
        self.assertEqual(len(first), 1)
        self.appendDp(w, 5)
        out = self.failedCycle(w, slot)
        self.assertEqual([e["name"] for e in w.spoolEntries()], first)
        self.assertIn("still unsent (will retry)", out)
        # The supervisor logs this one, so it has to name the condition.
        self.assertIn("still spooled", str(self.lastFailure))
        self.assertEqual(os.path.getsize(w.dpPath), 10 * RECORD_BYTES)

    def test_recovery_sends_the_spooled_delta_and_then_the_rest(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        _, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(10)))
        self.assertEqual(int(w.state["dpOffset"]), 10 * RECORD_BYTES)
        self.assertEqual(int(w.state["dpUploaded"]), 10)
        self.assertFalse(w.spoolPending())

    def test_records_stay_spooled_and_the_offset_does_not_move(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        self.assertEqual(int(w.state["dpOffset"]), 0)
        self.assertEqual(int(w.state.get("dpUploaded", 0)), 0)
        self.assertEqual(storedRecords(self.storeRoot), [])
        entries = w.spoolEntries()
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]["bytes"], 5 * RECORD_BYTES)
        self.assertEqual(entries[0]["key"], "dp/slot-%05d/%s-%016d-%s.bin"
                         % (slot, w.streamId, 0, sha256File(entries[0]["payload"])))
        self.assertEqual(os.path.basename(entries[0]["key"]), entries[0]["name"])
        with open(entries[0]["payload"], "rb") as fh:
            self.assertEqual(fh.read(), b"".join(record(i) for i in range(5)))
        self.assertTrue(w.spoolPending())

    def test_the_next_cycle_sends_the_spool_and_the_store_holds_each_record_once(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        self.appendDp(w, 3)
        out, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err, out)
        self.assertEqual(w.spoolEntries(), [])
        self.assertFalse(w.spoolPending())
        self.assertEqual(int(w.state["dpOffset"]), 8 * RECORD_BYTES)
        self.assertEqual(int(w.state["dpUploaded"]), 8)
        got = storedRecords(self.storeRoot)
        self.assertEqual(len(got), 8)
        self.assertEqual(sorted(got), sorted(record(i) for i in range(8)))
        self.assertEqual(len(set(got)), 8)

    def test_a_re_cut_of_the_same_records_lands_on_the_same_key(self):
        w, slot = self.claimed()
        self.appendDp(w, 4)
        self.failedCycle(w, slot)
        out, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err, out)
        self.assertEqual(len(storedKeys(self.storeRoot)), 1)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(4)))
        self.assertEqual(int(w.state["dpOffset"]), 4 * RECORD_BYTES)

    def test_a_restart_on_the_same_slot_credits_the_spooled_offset(self):
        """claim() prefers the slot just released; the leftover is that slot's.

        streamId is a new UUID in the new process, so crediting by streamId
        would leave dpOffset unmoved and the next cycle would recut [0, now)
        under a new key.
        """
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        w.slots.release(slot, w.owner)
        w.workLock.close()
        nxt = self.makeWorker()
        out, err = quiet(nxt.claimSlot)
        self.assertIsNone(err, out)
        self.assertEqual(int(nxt.state["slot"]), slot)
        self.assertNotEqual(nxt.streamId, w.streamId)
        self.assertFalse(nxt.spoolPending())
        self.assertEqual(int(nxt.state["dpOffset"]), 5 * RECORD_BYTES)
        self.appendDp(nxt, 3)
        out, err = quiet(nxt.uploadCycle, slot)
        self.assertIsNone(err, out)
        got = storedRecords(self.storeRoot)
        self.assertEqual(len(got), 8)
        self.assertEqual(sorted(got), sorted(record(i) for i in range(8)))
        self.assertEqual(len([k for k in storedKeys(self.storeRoot) if k.endswith(".bin")]), 2)

    def test_a_new_process_sends_what_the_old_one_left_under_the_old_slot(self):
        w, slot = self.claimed()
        self.appendDp(w, 6)
        self.failedCycle(w, slot)
        w.workLock.close()
        nxt = self.makeWorker()
        out, err = quiet(nxt.claimSlot)
        self.assertIsNone(err, out)
        self.assertNotEqual(int(nxt.state["slot"]), slot)
        self.assertEqual(nxt.spoolEntries(), [])
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(6)))
        for key in storedKeys(self.storeRoot):
            self.assertIn("slot-%05d" % slot, key)
        # The points belong to the slot they were cut for, not this one.
        self.assertEqual(int(nxt.state["dpOffset"]), 0)
        self.assertEqual(int(nxt.state.get("dpUploaded", 0)), 0)

    def test_reclaiming_the_original_slot_does_not_skip_the_new_file(self):
        """A leftover from a prior tenure must not move the new file's offset.

        claimSlot deletes dp.bin on a slot change.  If an outage then an
        intervening other-slot claim leave a same-slot spool entry on disk,
        crediting it by slot and offset=0 would jump past the start of the
        empty file the reclaimed slot is about to write.
        """
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        w.workLock.close()
        nxt = self.makeWorker()
        nxt.store = BrokenStore(nxt.store)
        out, err = quiet(nxt.claimSlot)
        self.assertIsNone(err, out)
        other = int(nxt.state["slot"])
        self.assertNotEqual(other, slot)
        self.assertTrue(nxt.spoolPending())
        self.assertEqual(int(nxt.state["dpOffset"]), 0)
        w.slots.release(slot, w.owner)
        nxt.slots.release(other, nxt.owner)
        nxt.workLock.close()
        nxt2 = self.makeWorker()
        out, err = quiet(nxt2.claimSlot)
        self.assertIsNone(err, out)
        self.assertEqual(int(nxt2.state["slot"]), slot)
        self.assertEqual(int(nxt2.state["dpOffset"]), 0)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))
        self.appendDp(nxt2, 3)
        out, err = quiet(nxt2.uploadCycle, slot)
        self.assertIsNone(err, out)
        self.assertEqual(int(nxt2.state["dpOffset"]), 3 * RECORD_BYTES)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(8)))
        self.assertEqual(len([k for k in storedKeys(self.storeRoot) if k.endswith(".bin")]), 2)


class AlreadyUploaded(SpoolCase):
    def test_a_key_the_store_already_has_is_dropped_without_re_uploading(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        # The put landed; only its acknowledgement was lost.
        dest = os.path.join(self.storeRoot, entry["key"])
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copyfile(entry["payload"], dest)
        counting = CountingStore(w.store)
        w.store = counting
        out, err = quiet(w.drainSpool, slot)
        self.assertIsNone(err, out)
        self.assertEqual(counting.puts, [])
        self.assertIn("already in the store", out)
        self.assertEqual(w.spoolEntries(), [])
        self.assertEqual(int(w.state["dpOffset"]), 5 * RECORD_BYTES)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))


class Budget(SpoolCase):
    def spoolThree(self, w, slot):
        """Three unsent deltas, as a box that has been through three processes.

        One process cannot stack three any more: uploadCycle refuses to cut a
        new delta while the spool holds one, so several entries mean several
        processes, or several slots, on the same disk.  They are spooled
        directly here for that reason.
        """
        entries = []
        for i in range(3):
            self.appendDp(w, 4)
            cut = os.path.join(w.work, "cut-%d.bin" % i)
            with open(cut, "wb") as fh:
                fh.write(b"".join(record(900 + 10 * i + j) for j in range(4)))
            w.streamId = "%032x" % (i + 1)
            key = "dp/slot-%05d/%s-%016d-%s.bin" % (
                slot, w.streamId, i * 4 * RECORD_BYTES, sha256File(cut))
            entries.append(w.spoolDelta(cut, key, slot, i * 4 * RECORD_BYTES))
            os.remove(cut)
        self.assertEqual(len(w.spoolEntries()), 3)
        return entries

    def test_the_oldest_entries_are_the_ones_kept(self):
        w, slot = self.claimed()
        entries = self.spoolThree(w, slot)
        w.cfg["spoolMaxBytes"] = w.entrySize(entries[0]) + w.entrySize(entries[1])
        out, err = quiet(w.enforceSpoolBudget)
        self.assertIsNone(err, out)
        kept = [e["name"] for e in w.spoolEntries()]
        self.assertEqual(kept, [entries[0]["name"], entries[1]["name"]])
        self.assertLessEqual(w.spoolBytes(), w.spoolBudget())

    def test_what_was_dropped_is_logged_with_its_record_count(self):
        w, slot = self.claimed()
        entries = self.spoolThree(w, slot)
        w.cfg["spoolMaxBytes"] = w.entrySize(entries[0])
        out, err = quiet(w.enforceSpoolBudget)
        self.assertIsNone(err, out)
        for gone in entries[1:]:
            self.assertIn(gone["name"], out)
            self.assertIn("%d records" % (gone["bytes"] // RECORD_BYTES), out)
            self.assertFalse(os.path.exists(gone["payload"]))
            self.assertFalse(os.path.exists(gone["manifest"]))
        self.assertIn("over budget", out)
        dropped = sum(e["bytes"] for e in entries[1:]) // RECORD_BYTES
        self.assertEqual(int(w.state["spoolDropped"]), dropped)

    def test_a_budget_that_fits_drops_nothing(self):
        w, slot = self.claimed()
        self.spoolThree(w, slot)
        w.cfg["spoolMaxBytes"] = w.spoolBytes()
        out, err = quiet(w.enforceSpoolBudget)
        self.assertIsNone(err, out)
        self.assertEqual(out, "")
        self.assertEqual(len(w.spoolEntries()), 3)
        self.assertNotIn("spoolDropped", w.state)


class Manifests(SpoolCase):
    """A strict campaign's dp envelope has to reach the store with its object."""

    def contract(self):
        return campaignContract({"storageProtocol": PROTOCOL, "curve": 131, "dpWeight": 34,
                                 "maxIters": 0, "packed": True, "workers": 385024, "batch": 16,
                                 "binarySha256": "a" * 64, "hostBinarySha256": "b" * 64,
                                 "sourceSha256": "c" * 64})

    def test_the_envelope_is_spooled_and_sent_with_its_delta(self):
        w, slot = self.claimed()
        w.contract = self.contract()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        self.assertTrue(entry["hasMeta"])
        self.assertTrue(os.path.exists(entry["payload"] + ".json"))
        out, err = quiet(w.drainSpool, slot)
        self.assertIsNone(err, out)
        self.assertFalse(w.spoolPending())
        blob = os.path.join(self.storeRoot, entry["key"])
        manifest = readJson(blob + ".json")
        verifyEnvelope(blob, manifest, w.contract, "dp")
        self.assertEqual(manifest["records"], 5)
        self.assertEqual(manifest["offset"], 0)


class Rotate(SpoolCase):
    def test_dp_file_is_kept_while_anything_is_unsent(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        # Pretend the offset is caught up: only the spool may hold rotation back.
        w.state["dpOffset"] = os.path.getsize(w.dpPath)
        w.saveState()
        out, err = quiet(w.rotateDpFile)
        self.assertIsNone(err, out)
        self.assertIn("not rotating dp.bin", out)
        self.assertTrue(os.path.exists(w.dpPath))

    def test_rotation_resumes_once_the_spool_is_empty(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        w.state["dpOffset"] = os.path.getsize(w.dpPath)
        w.saveState()
        out, err = quiet(w.drainSpool, slot)
        self.assertIsNone(err, out)
        self.assertFalse(w.spoolPending())
        stream = w.streamId
        _, err = quiet(w.rotateDpFile)
        self.assertIsNone(err)
        self.assertFalse(os.path.exists(w.dpPath))
        self.assertEqual(int(w.state["dpOffset"]), 0)
        self.assertNotEqual(w.streamId, stream)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))


class Fragments(SpoolCase):
    def test_a_payload_without_a_manifest_is_discarded_out_loud(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        os.remove(entry["manifest"])
        out, err = quiet(w.drainSpool, slot)
        self.assertIsNone(err, out)
        self.assertIn("discarding fragment", out)
        self.assertFalse(w.spoolPending())
        # dp.bin still holds them, so the next cut re-cuts them.
        self.assertEqual(int(w.state["dpOffset"]), 0)
        out, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err, out)
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))

    def test_a_manifest_without_a_payload_is_dropped_not_retried_forever(self):
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        os.remove(entry["payload"])
        counting = CountingStore(w.store)
        w.store = counting
        out, err = quiet(w.drainSpool, slot)
        self.assertIsNone(err, out)
        self.assertEqual(counting.puts, [])
        self.assertFalse(w.spoolPending())
        self.assertEqual(int(w.state["dpOffset"]), 0)

    def test_a_slot_change_sends_a_payload_that_never_got_a_manifest(self):
        """claimSlot deletes dp.bin; a fragment is then the only local copy."""
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        os.remove(entry["manifest"])
        self.assertTrue(os.path.exists(entry["payload"]))
        w.workLock.close()
        nxt = self.makeWorker()
        out, err = quiet(nxt.claimSlot)
        self.assertIsNone(err, out)
        self.assertNotEqual(int(nxt.state["slot"]), slot)
        self.assertEqual(nxt.spoolEntries(), [])
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))
        for key in storedKeys(self.storeRoot):
            self.assertIn("slot-%05d" % slot, key)
        self.assertEqual(int(nxt.state["dpOffset"]), 0)

    def test_retiring_sends_a_payload_that_never_got_a_manifest(self):
        """retireSlot clears state; the next claim must not sweep the fragment."""
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        os.remove(entry["manifest"])
        self.assertTrue(os.path.exists(entry["payload"]))
        _, err = quiet(w.retireSlot, slot, "checkpoint refused by the client")
        self.assertIsNone(err)
        self.assertEqual(w.state, {})
        w.workLock.close()
        nxt = self.makeWorker()
        out, err = quiet(nxt.claimSlot)
        self.assertIsNone(err, out)
        self.assertNotEqual(int(nxt.state["slot"]), slot)
        self.assertEqual(nxt.spoolEntries(), [])
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))
        for key in storedKeys(self.storeRoot):
            self.assertIn("slot-%05d" % slot, key)
        self.assertEqual(int(nxt.state["dpOffset"]), 0)

    def test_a_missing_state_file_still_sends_a_fragment(self):
        """state.json gone, checkpoint still names the slot the fragment belongs to."""
        w, slot = self.claimed()
        self.appendDp(w, 5)
        self.failedCycle(w, slot)
        entry = w.spoolEntries()[0]
        os.remove(entry["manifest"])
        self.assertTrue(os.path.exists(entry["payload"]))
        with open(w.ckptPath, "wb") as fh:
            fh.write(struct.pack("<8s6IQ", CKPT_MAGIC, 1, 131, 2, 4, 64, slot + 1, 1))
        os.remove(w.statePath)
        w.workLock.close()
        nxt = self.makeWorker()
        self.assertEqual(nxt.state, {})
        out, err = quiet(nxt.claimSlot)
        self.assertIsNone(err, out)
        self.assertNotEqual(int(nxt.state["slot"]), slot)
        self.assertEqual(nxt.spoolEntries(), [])
        self.assertEqual(sorted(storedRecords(self.storeRoot)),
                         sorted(record(i) for i in range(5)))
        for key in storedKeys(self.storeRoot):
            self.assertIn("slot-%05d" % slot, key)
        self.assertEqual(int(nxt.state["dpOffset"]), 0)


class V2Corpus(SpoolCase):
    """The witness corpus through the same spool, which v1 alone never tests.

    Every other case in this file writes a headerless 32-byte stream.  That
    left the whole format this client now writes by default uncovered, and it
    is where the offset accounting is hardest: a v2 cut starts at the file's
    16-byte header rather than at 0, so `dpOffset` of 0 and an entry at 16
    name the same place.  Comparing them for equality made the first upload of
    every v2 corpus fail to credit -- each cycle recut the whole file, dp.bin
    never rotated, and the merge ingested overlapping deltas.
    """

    def appendDpV2(self, w, count):
        fresh = not os.path.exists(w.dpPath) or os.path.getsize(w.dpPath) == 0
        with open(w.dpPath, "ab") as fh:
            if fresh:
                fh.write(v2Header())
            for i in range(self.written, self.written + count):
                fh.write(recordV2(i))
        self.written += count

    def test_the_first_cycle_credits_past_the_header(self):
        """The bug, stated as a test: dpOffset has to move off 0."""
        w, slot = self.claimed()
        self.appendDpV2(w, 5)
        out, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err, out)
        self.assertEqual(int(w.state["dpOffset"]),
                         DP_HEADER_BYTES + 5 * RECORD_BYTES_V2)
        self.assertEqual(int(w.state["dpUploaded"]), 5)
        self.assertEqual(w.spoolEntries(), [])
        self.assertEqual(sorted(storedRecordsV2(self.storeRoot)),
                         sorted(recordV2(i) for i in range(5)))

    def test_a_second_cycle_sends_only_what_is_new(self):
        """An offset that did not move resends the whole file every cycle."""
        w, slot = self.claimed()
        self.appendDpV2(w, 4)
        _, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err)
        firstKeys = storedKeys(self.storeRoot)
        self.appendDpV2(w, 3)
        out, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err, out)
        got = storedRecordsV2(self.storeRoot)
        self.assertEqual(len(got), 7, "a record was sent twice")
        self.assertEqual(len(set(got)), 7)
        self.assertEqual(sorted(got), sorted(recordV2(i) for i in range(7)))
        self.assertEqual(len(storedKeys(self.storeRoot)), len(firstKeys) + 1)
        self.assertEqual(int(w.state["dpOffset"]),
                         DP_HEADER_BYTES + 7 * RECORD_BYTES_V2)

    def test_every_delta_announces_the_format(self):
        """Deltas are standalone objects, so each carries its own header.

        Without this only a slot's first delta would announce v2 and every
        later one would be read as v1 -- invisible until the merge reports
        orbits nobody walked.
        """
        w, slot = self.claimed()
        for _ in range(3):
            self.appendDpV2(w, 2)
            _, err = quiet(w.uploadCycle, slot)
            self.assertIsNone(err)
        keys = storedKeys(self.storeRoot)
        self.assertEqual(len(keys), 3)
        for key in keys:
            with open(os.path.join(self.storeRoot, key), "rb") as fh:
                head = fh.read(DP_HEADER_BYTES)
            self.assertEqual(head[:len(DP_MAGIC_V2)], DP_MAGIC_V2, key)
            self.assertEqual(struct.unpack("<II", head[len(DP_MAGIC_V2):]),
                             (2, RECORD_BYTES_V2), key)
            size = os.path.getsize(os.path.join(self.storeRoot, key))
            self.assertEqual((size - DP_HEADER_BYTES) % RECORD_BYTES_V2, 0, key)

    def test_the_offset_advances_by_records_not_payload_bytes(self):
        """A delta's payload includes its header; the source stream does not.

        Advancing by payload bytes would step the offset past records that
        were never sent, once per cycle.
        """
        w, slot = self.claimed()
        self.appendDpV2(w, 6)
        _, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err)
        entries = storedKeys(self.storeRoot)
        self.assertEqual(len(entries), 1)
        payload = os.path.getsize(os.path.join(self.storeRoot, entries[0]))
        self.assertEqual(payload, DP_HEADER_BYTES + 6 * RECORD_BYTES_V2)
        # Off by exactly the header if payload bytes were used.
        self.assertEqual(int(w.state["dpOffset"]),
                         DP_HEADER_BYTES + 6 * RECORD_BYTES_V2)
        self.assertNotEqual(int(w.state["dpOffset"]),
                            DP_HEADER_BYTES + payload)

    def test_a_spooled_v2_delta_credits_on_recovery(self):
        """The failure path: the clamped offset has to reconcile there too."""
        w, slot = self.claimed()
        self.appendDpV2(w, 5)
        self.failedCycle(w, slot)
        self.assertEqual(int(w.state["dpOffset"]), 0)
        self.assertEqual(len(w.spoolEntries()), 1)
        self.assertEqual(w.spoolEntries()[0]["head"], DP_HEADER_BYTES)
        self.assertEqual(w.spoolEntries()[0]["stride"], RECORD_BYTES_V2)
        self.appendDpV2(w, 2)
        out, err = quiet(w.uploadCycle, slot)
        self.assertIsNone(err, out)
        self.assertEqual(w.spoolEntries(), [])
        got = storedRecordsV2(self.storeRoot)
        self.assertEqual(sorted(got), sorted(recordV2(i) for i in range(7)))
        self.assertEqual(len(set(got)), 7)
        self.assertEqual(int(w.state["dpOffset"]),
                         DP_HEADER_BYTES + 7 * RECORD_BYTES_V2)


if __name__ == "__main__":
    unittest.main()
