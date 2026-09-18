#!/usr/bin/env python3
"""Supervise one GPU's share of the ECC2K-130 campaign.

One process per GPU (systemd runs `worker@N` for every GPU on the instance).
It claims a *slot* -- a campaign-wide walk identity -- runs the client on that
slot, and keeps the slot's state durable so any instance can pick it up later:

  * seeds: the client derives every walk seed from --run-id, so a slot is a
    run id (slot + 1, run ids are 16-bit) and two live workers must never share
    one.  Slots are leased through a DynamoDB table; a lease outlives a crash
    by a few minutes and then the slot is free for the next claimant.
  * checkpoint: the client writes it locally on a timer; the supervisor copies
    it to S3 after every timer tick.  A replacement worker resumes exactly the
    walks in flight instead of throwing them away (about a quarter of all work
    at any instant sits in unreported walks).
  * distinguished points: the client appends 32-byte records to a local file;
    the supervisor uploads each new stretch of whole records as one immutable
    S3 object.  Points are uploaded *before* the checkpoint that follows them,
    so a resume can only re-report a point, never lose one.
  * a stretch that cannot be uploaded is kept in a local spool directory under
    its final object key and retried every cycle, so a store or network outage
    costs a delay rather than points.  The spool is a *process and transport*
    safety net only: it lives on the instance's disk, so a Spot reclamation
    takes anything still in it.  campaign.json's spoolMaxBytes caps it.
  * an exit code 6 (checkpoint refused) retires the slot rather than restarting
    it from scratch: re-walking a run id's seeds from their start points would
    repeat work already done and reported.

The client is stopped with SIGTERM, never killed, and finishes its launch,
flushes points and checkpoints before exiting.  Spot interruptions and
`systemctl stop` both arrive as SIGTERM here and are forwarded.

Environment (written to /etc/ecc2k130.env by bootstrap.sh):
  ECC_BUCKET, AWS_DEFAULT_REGION   the campaign bucket (slots live in it too)
  ECC_TABLE      optional DynamoDB table for slots instead of S3 objects
  ECC_GPU        GPU index on this instance (default 0)
  ECC_ROOT       directory holding campaign.json, the client, per-GPU work dirs
  ECC_CLIENT     client binary (default ECC_ROOT/ecc2k130)
  ECC_LOCAL_STORE  directory that stands in for S3 and DynamoDB (rehearsals)

No type hints, camelCase identifiers (project convention).
"""

import fcntl
import json
import os
import queue
import re
import shutil
import signal
import socket
import struct
import subprocess
import sys
import threading
import time
import urllib.request
import uuid

from protocol import (atomicJson, bindDirectory, campaignContract, envelope,
                      sha256File, verifyEnvelope)

# Frozen campaign fields a kernel rollout must not move. Duplicated from
# rollout.py so an already-booted box can pick up a new worker.py without
# also fetching that helper. kernelProtocol versions the client pointer
# and is not storageProtocol.
FROZEN_CAMPAIGN = (
    "curve", "dpWeight", "workers", "batch", "blockThreads", "minBlocks", "walk",
)
KERNEL_PROTOCOL = "ecc2k-kernel-v1"

PROGRESS_RE = re.compile(
    r"([\d.]+)\s+s\s+([\d.]+)\s+M it/s\s+(\d+)\s+iterations\s+(\d+)\s+dp\s+(\d+)\s+stored"
    r"(?:\s+(\d+)\s+dropped)?")
# The client's opening line names the grid it actually built.  Ada slots omit
# --threads and let autoThreads size it (usesCampaignWorkers), so a slot's walk
# count is not always campaign.json's workers x batch; the dashboard needs this
# number per slot to turn a checkpointed iteration base into group operations.
BANNER_RE = re.compile(r"=\s*(\d+)\s+walks,\s*dp weight")
RECORD_BYTES = 32
SPOOL_DIR = "spool"
SPOOL_MANIFEST = ".spool.json"
SPOOL_MAX_BYTES = 2 * 1024 * 1024 * 1024   # unsent points a worker may hold
CKPT_MAGIC = b"ECC2K130"
CKPT_ITER_OFFSET = 32      # magic[8] + version, m, threads, batch, lanes, runId (u32 each)
LEASE_SECONDS = 180
HEARTBEAT_SECONDS = 60
GRACE_SECONDS = 600        # the client checkpoints on SIGTERM; give it this long
MAX_SLOT = 65534           # run id = slot + 1 must fit in 16 bits


def log(msg):
    print(time.strftime("%Y-%m-%dT%H:%M:%SZ ", time.gmtime()) + msg, flush=True)


def sh(cmd, check=True):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError("%s failed (%d): %s" % (" ".join(cmd), r.returncode, (r.stderr or r.stdout).strip()))
    return r


def readJson(path, default=None):
    if not os.path.exists(path):
        return default
    with open(path) as fh:
        return json.load(fh)


def writeJson(path, obj):
    atomicJson(path, obj)


def copyFsync(src, dest):
    """Copy, and do not return until the bytes are on the disk.

    The spool's manifest is written with atomicJson, which fsyncs; without
    this the manifest could reach the disk first and claim records that were
    still only in the page cache.
    """
    with open(src, "rb") as fh, open(dest, "wb") as out:
        shutil.copyfileobj(fh, out)
        out.flush()
        os.fsync(out.fileno())


def checkpointIter(path):
    """iterBase from a checkpoint header, or -1 when the file is not one."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(CKPT_ITER_OFFSET + 8)
    except OSError:
        return -1
    if len(head) < CKPT_ITER_OFFSET + 8 or head[:8] != CKPT_MAGIC:
        return -1
    return struct.unpack_from("<Q", head, CKPT_ITER_OFFSET)[0]


def instanceId():
    try:
        req = urllib.request.Request("http://169.254.169.254/latest/api/token", method="PUT",
                                     headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"})
        token = urllib.request.urlopen(req, timeout=1).read().decode()
        req = urllib.request.Request("http://169.254.169.254/latest/meta-data/instance-id",
                                     headers={"X-aws-ec2-metadata-token": token})
        return urllib.request.urlopen(req, timeout=1).read().decode()
    except Exception:
        return socket.gethostname()


def gpuName(gpu):
    # Non-GPU clients (the FPGA host program) have no nvidia-smi; the
    # bootstrap tells us what the device is instead.
    if os.environ.get("ECC_DEVICE_NAME"):
        return os.environ["ECC_DEVICE_NAME"]
    try:
        r = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader", "-i", str(gpu)],
                           capture_output=True, text=True)
    except OSError:
        return "cpu"
    return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else "cpu"


# campaign.json workers=385024 is the RTX PRO 6000 / g7e preset. Ada (g6 L4,
# g6e L40S) auto-sizes; a checkpoint is only loadable into the same worker
# count, so Ada must not resume a Blackwell slot and g6 must not resume g6e.
BLACKWELL_FAMILIES = frozenset({"g7", "g7e"})
ADA_FAMILIES = frozenset({"g6", "g6e"})
LOCAL_FAMILIES = frozenset({"", "local", "cpu", None})


def instanceType():
    if os.environ.get("ECC_INSTANCE_TYPE"):
        return os.environ["ECC_INSTANCE_TYPE"]
    try:
        req = urllib.request.Request("http://169.254.169.254/latest/api/token", method="PUT",
                                     headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"})
        token = urllib.request.urlopen(req, timeout=1).read().decode()
        req = urllib.request.Request("http://169.254.169.254/latest/meta-data/instance-type",
                                     headers={"X-aws-ec2-metadata-token": token})
        return urllib.request.urlopen(req, timeout=1).read().decode()
    except Exception:
        return ""


def gpuFamily(name="", instance_type=""):
    """EC2 family used to pin slots and decide --threads vs autoThreads."""
    it = (instance_type or "").split(".")[0].lower()
    if it in ("g6", "g6e", "g7", "g7e", "g4dn", "g5", "g5g"):
        return it
    n = (name or "").lower()
    if "l40s" in n:
        return "g6e"
    if "rtx pro 6000" in n or "rtx 6000" in n:
        return "g7e"
    if "rtx pro 4500" in n or "rtx 4500" in n:
        return "g7"
    if "tesla t4" in n or n.endswith(" t4") or n == "t4":
        return "g4dn"
    if "l4" in n:
        return "g6"
    if n in ("cpu", "local"):
        return n
    return ""


def slotFamilyCompatible(slot_family, worker_family):
    """True if this worker may resume (or first-claim) the slot.

    Untagged slots are the live Blackwell corpus. Ada auto-sizes and would
    refuse those checkpoints (exit 6), retiring the run id, so Ada creates
    new slots instead. Rehearsal families stay compatible with anything.
    """
    if worker_family in LOCAL_FAMILIES:
        return True
    if not slot_family:
        return worker_family in BLACKWELL_FAMILIES
    return slot_family == worker_family


def usesCampaignWorkers(family):
    """Ada omits --threads so packedengine.autoThreads sizes the grid."""
    return family not in ADA_FAMILIES


def frozenCampaignMoved(current, nxt):
    """First frozen campaign field that changed, or None."""
    for key in FROZEN_CAMPAIGN:
        if key in current and key in nxt and current[key] != nxt[key]:
            return key
    return None


def campaignPointerMoved(current, nxt):
    """True when the bucket names a different client. Geometry moves are not a pointer move."""
    if frozenCampaignMoved(current, nxt):
        return False
    liveKey = current.get("binaryKey") or ""
    liveSha = current.get("binarySha256") or ""
    newKey = nxt.get("binaryKey") or ""
    newSha = nxt.get("binarySha256") or ""
    if newKey and newKey != liveKey:
        return True
    if newSha and newSha != liveSha:
        return True
    liveVer = int(current.get("kernelVersion") or 0)
    newVer = int(nxt.get("kernelVersion") or 0)
    if newVer and newVer != liveVer:
        return True
    if (nxt.get("kernelProtocol") or "") and (nxt.get("kernelProtocol") != (current.get("kernelProtocol") or "")):
        return True
    return False


# ---------------------------------------------------------------------------
# object store: S3, or a directory for rehearsals
# ---------------------------------------------------------------------------
class S3Store:
    def __init__(self, bucket):
        self.bucket = bucket

    def exists(self, key):
        r = subprocess.run(["aws", "s3api", "head-object", "--bucket", self.bucket, "--key", key],
                           capture_output=True, text=True)
        if r.returncode == 0:
            return True
        if "404" in r.stderr or "NoSuchKey" in r.stderr or "Not Found" in r.stderr:
            return False
        raise RuntimeError("object lookup failed; refusing to treat an access error as absence")

    def get(self, key, dest):
        if not self.exists(key):
            return False
        sh(["aws", "s3", "cp", "s3://%s/%s" % (self.bucket, key), dest, "--only-show-errors"])
        return True

    def put(self, src, key):
        sh(["aws", "s3", "cp", src, "s3://%s/%s" % (self.bucket, key), "--only-show-errors"])


class LocalStore:
    def __init__(self, root):
        self.root = root

    def exists(self, key):
        return os.path.exists(os.path.join(self.root, key))

    def get(self, key, dest):
        src = os.path.join(self.root, key)
        if not os.path.exists(src):
            return False
        shutil.copyfile(src, dest)
        return True

    def put(self, src, key):
        dest = os.path.join(self.root, key)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copyfile(src, dest + ".tmp")
        os.replace(dest + ".tmp", dest)


# ---------------------------------------------------------------------------
# slot leases: DynamoDB, or a locked JSON file for rehearsals
# ---------------------------------------------------------------------------
def dv(v):
    """DynamoDB attribute value for a Python scalar."""
    if isinstance(v, bool):
        return {"BOOL": v}
    if isinstance(v, (int, float)):
        return {"N": repr(v) if isinstance(v, float) else str(v)}
    return {"S": str(v)}


def fromDv(v):
    if "N" in v:
        n = v["N"]
        return float(n) if "." in n or "e" in n else int(n)
    if "S" in v:
        return v["S"]
    if "BOOL" in v:
        return v["BOOL"]
    return None


class DynamoSlots:
    def __init__(self, table):
        self.table = table

    def _run(self, *args):
        r = subprocess.run(["aws", "dynamodb"] + list(args) + ["--output", "json"],
                           capture_output=True, text=True)
        return r

    def scan(self):
        r = self._run("scan", "--table-name", self.table,
                      "--projection-expression", "#s, leaseUntil, #st, #o",
                      "--expression-attribute-names", json.dumps({"#s": "slot", "#st": "state", "#o": "owner"}))
        if r.returncode != 0:
            raise RuntimeError("dynamodb scan failed: " + r.stderr.strip())
        items = json.loads(r.stdout).get("Items", [])
        return [{k: fromDv(v) for k, v in it.items()} for it in items]

    def _update(self, slot, expr, names, values, condition=None):
        args = ["update-item", "--table-name", self.table, "--key", json.dumps({"slot": dv(slot)}),
                "--update-expression", expr, "--expression-attribute-names", json.dumps(names),
                "--expression-attribute-values", json.dumps(values)]
        if condition:
            args += ["--condition-expression", condition]
        r = self._run(*args)
        if r.returncode == 0:
            return True
        if "ConditionalCheckFailedException" in r.stderr:
            return False
        raise RuntimeError("dynamodb update failed: " + r.stderr.strip())

    def claim(self, owner, info):
        now = int(time.time())
        items = self.scan()
        free = sorted(it["slot"] for it in items
                      if it.get("state") not in ("retired", "solved", "error")
                      and int(it.get("leaseUntil") or 0) < now)
        names = {"#o": "owner", "#st": "state"}
        for slot in free:
            it = next((x for x in items if x.get("slot") == slot), {})
            if not slotFamilyCompatible(it.get("gpuFamily"), info.get("gpuFamily")):
                continue
            values = {":me": dv(owner), ":t": dv(now + LEASE_SECONDS), ":now": dv(now), ":active": dv("active"),
                      ":idle": dv("idle"), ":inst": dv(info["instance"]), ":gpu": dv(info["gpu"]),
                      ":gpuName": dv(info["gpuName"]), ":gpuFamily": dv(info.get("gpuFamily") or "")}
            # Only an expired or released lease may be taken, and only from a
            # slot that is still walking: retired, solved and error slots keep
            # their run id forever so its seeds are never walked twice.
            ok = self._update(slot, "SET #o = :me, leaseUntil = :t, claimedAt = :now, #st = :active, "
                              "instance = :inst, gpu = :gpu, gpuName = :gpuName, gpuFamily = :gpuFamily", names, values,
                              "(attribute_not_exists(leaseUntil) OR leaseUntil < :now) AND "
                              "(attribute_not_exists(#st) OR #st = :active OR #st = :idle)")
            if ok:
                return slot
        nextSlot = (max((it["slot"] for it in items), default=-1)) + 1
        for _ in range(64):
            if nextSlot > MAX_SLOT:
                raise RuntimeError("run ids exhausted")
            item = {"slot": dv(nextSlot), "owner": dv(owner), "leaseUntil": dv(now + LEASE_SECONDS),
                    "claimedAt": dv(now), "createdAt": dv(now), "state": dv("active"),
                    "instance": dv(info["instance"]), "gpu": dv(info["gpu"]),
                    "gpuName": dv(info["gpuName"]), "gpuFamily": dv(info.get("gpuFamily") or "")}
            r = self._run("put-item", "--table-name", self.table, "--item", json.dumps(item),
                          "--condition-expression", "attribute_not_exists(slot)")
            if r.returncode == 0:
                return nextSlot
            if "ConditionalCheckFailedException" not in r.stderr:
                raise RuntimeError("dynamodb put failed: " + r.stderr.strip())
            nextSlot += 1
        raise RuntimeError("could not allocate a slot")

    def heartbeat(self, slot, owner, fields):
        now = int(time.time())
        expr = ["#o = :me", "leaseUntil = :t", "updatedAt = :now"]
        names = {"#o": "owner"}
        values = {":me": dv(owner), ":t": dv(now + LEASE_SECONDS), ":now": dv(now)}
        for i, (k, v) in enumerate(sorted(fields.items())):
            names["#f%d" % i] = k
            values[":v%d" % i] = dv(v)
            expr.append("#f%d = :v%d" % (i, i))
        return self._update(slot, "SET " + ", ".join(expr), names, values,
                            "#o = :me AND leaseUntil >= :now")

    def release(self, slot, owner, state="idle", extra=None):
        names = {"#o": "owner", "#st": "state"}
        values = {":me": dv(owner), ":zero": dv(0), ":st": dv(state)}
        expr = ["leaseUntil = :zero", "#st = :st"]
        for i, (k, v) in enumerate(sorted((extra or {}).items())):
            names["#f%d" % i] = k
            values[":v%d" % i] = dv(v)
            expr.append("#f%d = :v%d" % (i, i))
        return self._update(slot, "SET " + ", ".join(expr), names, values, "#o = :me")


class S3Slots:
    """Same contract as DynamoSlots on S3 objects with conditional writes.

    One object per slot, slots/slot-NNNNN.json.  A claim, heartbeat or release
    is a read followed by a PutObject with If-Match on the ETag just read, so
    two workers racing for one slot cannot both win; a new slot is a PutObject
    with If-None-Match: *.  S3 reads are strongly consistent, which is what
    makes the read-modify-write safe.  Needs only s3:GetObject/PutObject/
    ListBucket, so the instance role stays minimal."""

    def __init__(self, bucket):
        self.bucket = bucket
        self.tmp = os.path.join(os.environ.get("TMPDIR", "/tmp"), "ecc-slot-%d.json" % os.getpid())

    def _key(self, slot):
        return "slots/slot-%05d.json" % slot

    def _get(self, slot):
        r = subprocess.run(["aws", "s3api", "get-object", "--bucket", self.bucket, "--key", self._key(slot),
                            self.tmp, "--output", "json"], capture_output=True, text=True)
        if r.returncode != 0:
            if "NoSuchKey" in r.stderr or "404" in r.stderr or "Not Found" in r.stderr:
                return None, None
            raise RuntimeError("s3 get-object failed: " + r.stderr.strip())
        etag = json.loads(r.stdout)["ETag"]
        with open(self.tmp) as fh:
            return json.load(fh), etag

    def _put(self, slot, item, ifMatch=None, create=False):
        writeJson(self.tmp, item)
        cmd = ["aws", "s3api", "put-object", "--bucket", self.bucket, "--key", self._key(slot),
               "--body", self.tmp, "--content-type", "application/json"]
        if create:
            cmd += ["--if-none-match", "*"]
        elif ifMatch:
            cmd += ["--if-match", ifMatch]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            return True
        if "PreconditionFailed" in r.stderr or "412" in r.stderr or "ConditionalRequestConflict" in r.stderr:
            return False
        raise RuntimeError("s3 put-object failed: " + r.stderr.strip())

    def scan(self):
        r = subprocess.run(["aws", "s3api", "list-objects-v2", "--bucket", self.bucket, "--prefix", "slots/",
                            "--query", "Contents[].Key", "--output", "json"], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError("s3 list failed: " + r.stderr.strip())
        keys = json.loads(r.stdout) or []
        items = []
        for key in keys:
            name = os.path.basename(key)
            if not (name.startswith("slot-") and name.endswith(".json")):
                continue
            slot = int(name[5:-5])
            item, etag = self._get(slot)
            if item is not None:
                item["slot"] = slot
                item["_etag"] = etag
                items.append(item)
        return items

    def claim(self, owner, info):
        now = int(time.time())
        items = self.scan()
        for it in sorted(items, key=lambda x: x["slot"]):
            if it.get("state") not in (None, "active", "idle") or int(it.get("leaseUntil") or 0) >= now:
                continue
            if not slotFamilyCompatible(it.get("gpuFamily"), info.get("gpuFamily")):
                continue
            etag, slot = it["_etag"], it["slot"]
            candidate = {k: v for k, v in it.items() if k not in ("_etag", "slot")}
            candidate.update(owner=owner, leaseUntil=now + LEASE_SECONDS, claimedAt=now, state="active", **info)
            if self._put(slot, candidate, ifMatch=etag):
                return slot
        nextSlot = max((it["slot"] for it in items), default=-1) + 1
        for _ in range(64):
            if nextSlot > MAX_SLOT:
                raise RuntimeError("run ids exhausted")
            item = dict(owner=owner, leaseUntil=now + LEASE_SECONDS, claimedAt=now, createdAt=now,
                        state="active", **info)
            if self._put(nextSlot, item, create=True):
                return nextSlot
            nextSlot += 1
        raise RuntimeError("could not allocate a slot")

    def _modify(self, slot, owner, fn):
        item, etag = self._get(slot)
        if item is None or item.get("owner") != owner or item.get("leaseUntil", 0) < int(time.time()):
            return False
        fn(item)
        return self._put(slot, item, ifMatch=etag)

    def heartbeat(self, slot, owner, fields):
        now = int(time.time())
        return self._modify(slot, owner, lambda it: it.update(fields, leaseUntil=now + LEASE_SECONDS, updatedAt=now))

    def release(self, slot, owner, state="idle", extra=None):
        return self._modify(slot, owner, lambda it: it.update(extra or {}, leaseUntil=0, state=state))


class LocalSlots:
    """Same contract as DynamoSlots on one JSON file under an exclusive lock."""

    def __init__(self, path):
        self.path = path

    def _locked(self, fn):
        with open(self.path + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            items = readJson(self.path, {})
            out = fn(items)
            writeJson(self.path, items)
            return out

    def scan(self):
        return [dict(it, slot=int(k)) for k, it in readJson(self.path, {}).items()]

    def claim(self, owner, info):
        def fn(items):
            now = int(time.time())
            for k in sorted(items, key=int):
                it = items[k]
                if it.get("state") in ("retired", "solved", "error") or int(it.get("leaseUntil") or 0) >= now:
                    continue
                if not slotFamilyCompatible(it.get("gpuFamily"), info.get("gpuFamily")):
                    continue
                it.update(owner=owner, leaseUntil=now + LEASE_SECONDS, claimedAt=now, state="active", **info)
                return int(k)
            slot = max((int(k) for k in items), default=-1) + 1
            items[str(slot)] = dict(owner=owner, leaseUntil=now + LEASE_SECONDS, claimedAt=now,
                                    createdAt=now, state="active", **info)
            return slot
        return self._locked(fn)

    def heartbeat(self, slot, owner, fields):
        def fn(items):
            it = items.get(str(slot))
            if not it or it.get("owner") != owner or it.get("leaseUntil", 0) < int(time.time()):
                return False
            it.update(fields, leaseUntil=int(time.time()) + LEASE_SECONDS, updatedAt=int(time.time()))
            return True
        return self._locked(fn)

    def release(self, slot, owner, state="idle", extra=None):
        def fn(items):
            it = items.get(str(slot))
            if not it or it.get("owner") != owner:
                return False
            it.update(extra or {}, leaseUntil=0, state=state)
            return True
        return self._locked(fn)


# ---------------------------------------------------------------------------
class Worker:
    def __init__(self):
        self.root = os.environ.get("ECC_ROOT", os.getcwd())
        self.gpu = int(os.environ.get("ECC_GPU", "0"))
        self.client = os.environ.get("ECC_CLIENT", os.path.join(self.root, "ecc2k130"))
        local = os.environ.get("ECC_LOCAL_STORE")
        if local:
            os.makedirs(local, exist_ok=True)
            self.store = LocalStore(local)
            self.slots = LocalSlots(os.path.join(local, "slots.json"))
        else:
            self.store = S3Store(os.environ["ECC_BUCKET"])
            # DynamoDB when a table is named, otherwise slots live in the
            # bucket itself (fewer services, fewer permissions).
            table = os.environ.get("ECC_TABLE", "")
            self.slots = DynamoSlots(table) if table else S3Slots(os.environ["ECC_BUCKET"])
        self.instance = instanceId()
        self.owner = "%s:gpu%d:%s" % (self.instance, self.gpu, uuid.uuid4().hex)
        self.work = os.path.join(self.root, "gpu%d" % self.gpu)
        os.makedirs(self.work, exist_ok=True)
        self.statePath = os.path.join(self.work, "state.json")
        self.state = readJson(self.statePath, {})
        self.stopping = False
        self.proc = None
        self.cfg = None
        self.contract = None
        self.leaseLost = False
        self.lastBeatSuccess = None
        self.streamId = uuid.uuid4().hex
        # Prevent two local supervisors sharing offsets, checkpoint paths or a GPU.
        self.workLock = open(os.path.join(self.work, "worker.lock"), "a+")
        fcntl.flock(self.workLock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.gpuName = gpuName(self.gpu)
        self.instanceType = instanceType()
        self.gpuFamily = gpuFamily(self.gpuName, self.instanceType)
        signal.signal(signal.SIGTERM, self.onStop)
        signal.signal(signal.SIGINT, self.onStop)

    def onStop(self, signo, frame):
        if not self.stopping:
            log("signal %d: stopping gracefully" % signo)
        self.stopping = True
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.send_signal(signal.SIGTERM)
            except OSError:
                pass

    def saveState(self):
        writeJson(self.statePath, self.state)

    # ---- campaign configuration ------------------------------------------
    def loadConfig(self, verifyBinary=True):
        path = os.path.join(self.work, "campaign.json")
        if not self.store.get("campaign.json", path):
            raise RuntimeError("campaign.json missing from the store")
        self.cfg = readJson(path)
        for key in ("curve", "steps", "checkpointEvery"):
            if key not in self.cfg:
                raise RuntimeError("campaign.json lacks %r" % key)
        if self.cfg.get("storageProtocol"):
            self.contract = campaignContract(self.cfg)
            if verifyBinary and sha256File(self.client) != self.cfg["binarySha256"]:
                raise RuntimeError("client binary hash differs from campaign")
            bindDirectory(self.work, self.contract)
        elif not os.environ.get("ECC_ALLOW_LEGACY_STORAGE"):
            raise RuntimeError("unversioned campaign: set storageProtocol; legacy storage requires ECC_ALLOW_LEGACY_STORAGE=1")
        self.verifyKernelPin(verifyBinary)

    def verifyKernelPin(self, verifyBinary=True):
        """When kernelProtocol is live, the client must match the pinned hash.

        Independent of storageProtocol: the live corpus can adopt v1 without
        migrating points. Version 0 is an empty pointer (no binary yet).
        """
        proto = self.cfg.get("kernelProtocol")
        if not proto:
            return
        if proto != KERNEL_PROTOCOL:
            raise RuntimeError("unknown kernelProtocol %r" % proto)
        if int(self.cfg.get("kernelVersion") or 0) < 1 or not verifyBinary:
            return
        if self.gpuFamily in ("cpu", "local") or not self.cfg.get("packed", True):
            want = self.cfg.get("hostBinarySha256") or self.cfg.get("binarySha256")
        else:
            want = self.cfg.get("binarySha256")
        if not want:
            raise RuntimeError("kernelVersion %s has no pinned client hash" % self.cfg.get("kernelVersion"))
        if sha256File(self.client) != want:
            raise RuntimeError("client binary hash differs from kernelVersion %s" % self.cfg.get("kernelVersion"))

    def clientStoreKey(self):
        """S3 key of the executable this worker should run.

        GPU walkers take binaryKey. CPU walkers take hostBinaryKey so a
        kernel rollout does not hand them the CUDA client. Rehearsals use
        binaryKey=local and ECC_CLIENT; those do not re-fetch.
        """
        key = self.cfg.get("binaryKey") or ""
        if key in ("", "local"):
            return ""
        if self.gpuFamily in ("cpu", "local") or not self.cfg.get("packed", True):
            return self.cfg.get("hostBinaryKey") or key
        return key

    def fetchClient(self):
        key = self.clientStoreKey()
        if not key:
            return False
        tmp = self.client + ".next"
        if not self.store.get(key, tmp):
            raise RuntimeError("failed to download %s" % key)
        if key == (self.cfg.get("hostBinaryKey") or ""):
            want = self.cfg.get("hostBinarySha256") or ""
        else:
            want = self.cfg.get("binarySha256") or ""
        if want and sha256File(tmp) != want:
            os.remove(tmp)
            raise RuntimeError("downloaded client hash differs from campaign")
        os.chmod(tmp, 0o755)
        os.replace(tmp, self.client)
        return True

    def campaignPointerChanged(self):
        """True when campaign.json names a new client we are allowed to load."""
        remote = os.path.join(self.work, "campaign.json.next")
        try:
            if not self.store.get("campaign.json", remote):
                return False
            nxt = readJson(remote)
        except Exception as e:
            log("campaign pointer poll failed: %s" % e)
            return False
        moved = frozenCampaignMoved(self.cfg, nxt)
        if moved:
            log("refusing campaign.json: frozen field %s moved; keep walking the current client" % moved)
            return False
        if campaignPointerMoved(self.cfg, nxt):
            log("campaign client moved to %s (kernelVersion %s)"
                % (nxt.get("binaryKey") or nxt.get("binarySha256"), nxt.get("kernelVersion")))
            return True
        return False

    def reloadClient(self):
        """Re-fetch campaign.json and the client it names. Slot is kept."""
        log("reloading campaign.json and client")
        remote = os.path.join(self.work, "campaign.json.next")
        if not self.store.get("campaign.json", remote):
            raise RuntimeError("campaign.json missing during reload")
        nxt = readJson(remote)
        moved = frozenCampaignMoved(self.cfg, nxt)
        if moved:
            raise RuntimeError("frozen field %s moved; refusing reload" % moved)
        prev = self.cfg
        self.cfg = nxt
        try:
            self.fetchClient()
        except Exception:
            self.cfg = prev
            raise
        path = os.path.join(self.work, "campaign.json")
        os.replace(remote, path)
        self.loadConfig(verifyBinary=True)

    def clientCommand(self, slot):
        c = self.cfg
        cmd = [self.client, "--curve", str(c["curve"]), "--steps", str(c["steps"]), "--launches", "0",
               "--run-id", str(slot + 1), "--dp-file", self.dpPath, "--checkpoint", self.ckptPath,
               "--checkpoint-every", str(int(c["checkpointEvery"])), "--verify", str(int(c.get("verify", 0)))]
        if c.get("packed", False):
            cmd += ["--packed", "--device", "0"]
        if c.get("workers") and usesCampaignWorkers(self.gpuFamily):
            cmd += ["--threads", str(int(c["workers"]))]
        if c.get("dpWeight", -1) >= 0:
            cmd += ["--dp-weight", str(int(c["dpWeight"]))]
        if c.get("maxIters"):
            cmd += ["--max-iters", str(int(c["maxIters"]))]
        if c.get("loadMax"):
            cmd += ["--load-max", str(int(c["loadMax"]))]
        if c.get("dpCap"):
            cmd += ["--dp-cap", str(int(c["dpCap"]))]
        cmd += [str(a) for a in c.get("extraArgs", [])]
        return cmd

    # ---- slot lifecycle ---------------------------------------------------
    @property
    def dpPath(self):
        return os.path.join(self.work, "dp.bin")

    @property
    def ckptPath(self):
        return os.path.join(self.work, "walk.ck")

    def ckptKey(self, slot):
        return "ckpt/slot-%05d.ck" % slot

    def claimSlot(self):
        info = {"instance": self.instance, "gpu": self.gpu, "gpuName": self.gpuName,
                "gpuFamily": self.gpuFamily}
        slot = self.slots.claim(self.owner, info)
        self.lastBeatSuccess = time.monotonic()
        if self.contract:
            records = [r for r in self.slots.scan() if r["slot"] == slot]
            record = records[0]
            if record.get("campaignId") not in (None, self.contract["id"]):
                raise RuntimeError("slot belongs to a different campaign")
            if record.get("campaignId") is None and (record.get("ckptIter", -1) >= 0 or self.store.exists(self.ckptKey(slot))):
                raise RuntimeError("legacy slot checkpoint requires audited migration")
            if not self.slots.heartbeat(slot, self.owner, {"campaignId": self.contract["id"]}):
                raise RuntimeError("lease lost during campaign binding")
        log("claimed slot %d (run id %d)" % (slot, slot + 1))
        if self.state.get("slot") != slot:
            # A different slot than this work dir last held: the live dp.bin
            # and checkpoint do not apply.  Unsent spool entries do -- they
            # were cut for the slot they name.  Finish any payload-only
            # fragment and forget leftover offsets before deleting dp.bin,
            # so sweepSpool cannot drop the only copy and creditSpoolEntry
            # cannot treat those records as a prefix of the new empty file.
            old = self.state.get("slot")
            if old is not None:
                self.completeSpoolFragments(old)
            self.detachSpoolFromLiveFile()
            for name in ("dp.bin", "walk.ck", "walk.ck.snap", "walk.ck.remote"):
                p = os.path.join(self.work, name)
                if os.path.exists(p):
                    os.remove(p)
            self.state = {"slot": slot, "dpOffset": 0, "ckptIter": -1, "dpUploaded": 0, "walks": 0}
            self.saveState()
        # Resume from whichever checkpoint is further along: the one left here by
        # a previous run on this instance, or the one another instance uploaded.
        remote = self.ckptPath + ".remote"
        if self.contract:
            # The registry's conditional write is the checkpoint commit point.
            # A stale writer may leave immutable blobs but cannot change this pointer.
            key = record.get("checkpointKey")
            if key:
                meta = remote + ".json"
                if not self.store.get(key, remote) or not self.store.get(key + ".json", meta):
                    raise RuntimeError("committed checkpoint or manifest is missing")
                verifyEnvelope(remote, readJson(meta), self.contract, "checkpoint")
                os.replace(remote, self.ckptPath)
                log("downloaded checkpoint at iteration %d" % checkpointIter(self.ckptPath))
            # Do not trust an uncommitted local checkpoint after process restart.
            elif os.path.exists(self.ckptPath):
                raise RuntimeError("uncommitted local checkpoint requires recovery review")
        elif self.store.get(self.ckptKey(slot), remote):
            if checkpointIter(remote) > checkpointIter(self.ckptPath):
                os.replace(remote, self.ckptPath)
                log("downloaded checkpoint at iteration %d" % checkpointIter(self.ckptPath))
            else:
                os.remove(remote)
        if os.path.exists(self.ckptPath):
            log("local checkpoint at iteration %d" % checkpointIter(self.ckptPath))
        # Points a previous process on this disk cut but never sent are still
        # owed to the campaign, under the slot they were cut for.
        self.drainSpool(slot)
        return slot

    def retireSlot(self, slot, reason):
        log("retiring slot %d: %s" % (slot, reason))
        if os.path.exists(self.ckptPath):
            self.store.put(self.ckptPath, "ckpt/retired/slot-%05d.ck" % slot)
        self.slots.release(slot, self.owner, state="retired", extra={"reason": reason})
        self.state = {}
        self.saveState()

    # ---- local spool of unsent points -------------------------------------
    @property
    def spoolDir(self):
        return os.path.join(self.work, SPOOL_DIR)

    def spoolBudget(self):
        return int((getattr(self, "cfg", None) or {}).get("spoolMaxBytes", SPOOL_MAX_BYTES))

    def spoolEntries(self):
        """Cut deltas the store has not acknowledged, oldest first.

        Each entry keeps its own object key rather than deriving one: a spool
        left behind by a previous process holds points cut for the slot that
        process had, which is not necessarily the slot this one claimed, and
        re-filing them under the wrong slot would misattribute the walk.
        """
        if not os.path.isdir(self.spoolDir):
            return []
        out = []
        for name in sorted(os.listdir(self.spoolDir)):
            if not name.endswith(SPOOL_MANIFEST):
                continue
            entry = readJson(os.path.join(self.spoolDir, name))
            if not entry or not entry.get("key") or not entry.get("name"):
                continue
            entry["manifest"] = os.path.join(self.spoolDir, name)
            entry["payload"] = os.path.join(self.spoolDir, entry["name"])
            out.append(entry)
        out.sort(key=lambda e: (e.get("createdAt", 0), int(e.get("offset", 0)), e["name"]))
        return out

    def spoolBytes(self):
        total = 0
        if os.path.isdir(self.spoolDir):
            for name in os.listdir(self.spoolDir):
                path = os.path.join(self.spoolDir, name)
                if os.path.isfile(path):
                    total += os.path.getsize(path)
        return total

    def spoolPending(self):
        """True while anything at all sits in the spool, manifest or not."""
        return os.path.isdir(self.spoolDir) and bool(os.listdir(self.spoolDir))

    def entrySize(self, entry):
        total = 0
        for path in (entry["payload"], entry["payload"] + ".json", entry["manifest"]):
            if os.path.isfile(path):
                total += os.path.getsize(path)
        return total

    def spoolDelta(self, src, key, slot, offset, metaPath=None):
        """Hold a cut delta on disk under its final key until the store has it.

        The five-hour stall of 2026-09-17 was an ingest fault rather than an
        upload one, but it showed what an outage of the upload path would
        have cost: the only copy of an unsent stretch was dp.bin, which
        claimSlot deletes when the work dir changes slot and rotateDpFile
        deletes on a rollout restart.  The manifest is written last, so a
        payload without one is a fragment and never a record of work.
        """
        os.makedirs(self.spoolDir, exist_ok=True)
        name = os.path.basename(key)
        payload = os.path.join(self.spoolDir, name)
        copyFsync(src, payload + ".part")
        os.replace(payload + ".part", payload)
        if metaPath:
            copyFsync(metaPath, payload + ".json")
        entry = {"name": name, "key": key, "slot": slot, "streamId": self.streamId,
                 "offset": int(offset), "bytes": os.path.getsize(payload),
                 "hasMeta": bool(metaPath), "createdAt": time.time()}
        writeJson(payload + SPOOL_MANIFEST, entry)
        entry["manifest"] = payload + SPOOL_MANIFEST
        entry["payload"] = payload
        return entry

    def removeSpoolEntry(self, entry):
        # Manifest last: a crash here leaves an entry that drains as a missing
        # payload, never one that claims records the disk no longer holds.
        for path in (entry["payload"] + ".json", entry["payload"], entry["manifest"]):
            if os.path.exists(path):
                os.remove(path)

    def completeSpoolFragments(self, slot):
        """Give payload-only spool files a manifest for this slot.

        spoolDelta writes the payload first and the manifest last.  sweepSpool
        would drop a lone payload as a fragment whose records are still in
        dp.bin -- true only until claimSlot deletes that file on a slot
        change.  The payload's name is the basename of its final key.
        """
        if slot is None or not os.path.isdir(self.spoolDir):
            return
        claimed = {e["name"] for e in self.spoolEntries()}
        fallbackOffset = int(self.state.get("dpOffset", 0))
        for name in sorted(os.listdir(self.spoolDir)):
            if not name.endswith(".bin") or name in claimed:
                continue
            payload = os.path.join(self.spoolDir, name)
            if not os.path.isfile(payload):
                continue
            size = os.path.getsize(payload)
            if size <= 0 or size % RECORD_BYTES:
                continue
            parts = name[:-4].rsplit("-", 2)
            if len(parts) == 3 and parts[1].isdigit():
                streamId, offset = parts[0], int(parts[1])
            else:
                streamId, offset = self.streamId, fallbackOffset
            meta = payload + ".json"
            entry = {"name": name, "key": "dp/slot-%05d/%s" % (slot, name),
                     "slot": slot, "streamId": streamId, "offset": offset,
                     "bytes": size, "hasMeta": os.path.isfile(meta),
                     "createdAt": time.time()}
            writeJson(payload + SPOOL_MANIFEST, entry)

    def detachSpoolFromLiveFile(self):
        """Stop leftover entries from indexing a dp.bin that is about to go.

        creditSpoolEntry treats a matching slot and offset as the records
        still being in the live file.  After a slot change that file is a
        new empty one, and a later reclaim of the same slot number would
        jump dpOffset past its prefix.  Clearing the offset (not the
        payload) keeps the records publishable on their own key.
        """
        for entry in self.spoolEntries():
            body = {k: v for k, v in entry.items() if k not in ("manifest", "payload")}
            body["offset"] = -1
            writeJson(entry["manifest"], body)

    def sweepSpool(self):
        """Discard spool files that no manifest claims.

        A payload is written before dpOffset moves, so an interrupted write
        costs nothing: those records are still in dp.bin and the next cycle
        cuts them again.  rotateDpFile refusing to run while the spool is
        non-empty is what keeps that true on a rollout.  claimSlot changing
        slot completes any payload-only fragment before it deletes dp.bin,
        so a crash in spoolDelta cannot be swept away with the file.
        """
        if not os.path.isdir(self.spoolDir):
            return
        known = set()
        for entry in self.spoolEntries():
            known.update((entry["name"], entry["name"] + ".json",
                          entry["name"] + SPOOL_MANIFEST))
        for name in sorted(os.listdir(self.spoolDir)):
            path = os.path.join(self.spoolDir, name)
            if name in known or not os.path.isfile(path):
                continue
            log("spool: discarding fragment %s (%d bytes); its records are still in dp.bin"
                % (name, os.path.getsize(path)))
            os.remove(path)

    def creditSpoolEntry(self, entry, slot):
        """Move dpOffset past a delta the store has, and only then.

        The offset indexes one dp file, so only an entry for this slot whose
        offset is where this process is reading may move it.  streamId is not
        the test: it is a fresh UUID each process start and is not in
        state.json, so a same-slot restart would upload the leftover and then
        recut that prefix.  Rotation will not drop dp.bin while the spool
        holds anything, so a matching slot and offset are the file still here.
        claimSlot will drop it on a slot change, and forgets leftover offsets
        first so a later claim of the same slot number cannot treat those
        records as a prefix of the new empty file.  An entry for another slot
        is published on its own key and the offset stays where it is, because
        skipping the prefix of a file those records are not in would lose the
        points at the front of it.
        """
        if slot is None or entry.get("slot") != slot:
            return
        if int(entry.get("offset", -1)) != int(self.state.get("dpOffset", 0)):
            return
        self.state["dpOffset"] = int(entry["offset"]) + int(entry["bytes"])
        # Cumulative across dp file rotations, so the dashboard's count
        # is this slot's whole contribution.
        self.state["dpUploaded"] = int(self.state.get("dpUploaded", 0)) + int(entry["bytes"]) // RECORD_BYTES
        self.saveState()

    def uploadSpoolEntry(self, entry, slot=None, source=None, metaSource=None, check=True):
        """Put one delta and unspool it only once the store holds it.

        A key names the bytes it carries, so an attempt that failed after the
        object landed leaves the same key in place; treating that as sent
        keeps a retry from re-uploading a stretch the store already has.
        """
        payload = entry["payload"]
        src = payload if os.path.exists(payload) else source
        if src is None:
            log("spool: %s claims records that are not on disk; dropping the manifest" % entry["name"])
            if os.path.exists(entry["manifest"]):
                os.remove(entry["manifest"])
            return False
        meta = payload + ".json"
        metaSrc = meta if os.path.exists(meta) else (metaSource if entry.get("hasMeta") else None)
        if check and self.store.exists(entry["key"]):
            log("spool: %s is already in the store; dropping the local copy" % entry["name"])
            if metaSrc and not self.store.exists(entry["key"] + ".json"):
                self.store.put(metaSrc, entry["key"] + ".json")
        else:
            self.store.put(src, entry["key"])
            if metaSrc:
                self.store.put(metaSrc, entry["key"] + ".json")
        self.removeSpoolEntry(entry)
        self.creditSpoolEntry(entry, slot)
        return True

    def drainSpool(self, slot=None):
        """Send what earlier cycles, or an earlier process, could not.

        Never fatal: an unreachable store is the condition the spool exists
        for, so a failure is a log line and another attempt next cycle.  The
        first failure stops the drain because the cause is almost always the
        store itself, and the rest of the queue would only hammer it.
        """
        self.sweepSpool()
        entries = self.spoolEntries()
        if not entries:
            return 0
        records = sum(int(e.get("bytes", 0)) for e in entries) // RECORD_BYTES
        log("spool: draining %d unsent delta(s), %d records, %d bytes"
            % (len(entries), records, self.spoolBytes()))
        sent = 0
        for entry in entries:
            try:
                if self.uploadSpoolEntry(entry, slot):
                    sent += 1
            except Exception as e:
                log("spool: %s still unsent (will retry): %s" % (entry["name"], e))
                break
        return sent

    def enforceSpoolBudget(self):
        """Hold the spool to spoolMaxBytes by dropping the newest entries.

        A worker that cannot reach the store must not fill the disk out from
        under the client.  The oldest entries are the ones kept: they are the
        ones whose dp.bin may already be gone, while the newest were cut from
        the live file and the next cycle cuts them again.  What goes is said
        out loud, with its record count, and counted into spoolDropped, which
        the heartbeat publishes: a worker shedding points must not be a thing
        only its own log knows.
        """
        budget = self.spoolBudget()
        entries = self.spoolEntries()
        total = self.spoolBytes()
        dropped = 0
        while total > budget and entries:
            entry = entries.pop()
            count = int(entry.get("bytes", 0)) // RECORD_BYTES
            size = self.entrySize(entry)
            self.removeSpoolEntry(entry)
            log("spool over budget (%d > %d bytes): dropped %s, %d records no longer held locally"
                % (total, budget, entry["name"], count))
            total -= size
            dropped += count
        if dropped:
            self.state["spoolDropped"] = int(self.state.get("spoolDropped", 0)) + dropped
            self.saveState()
            log("spool: %d records dropped from the spool on this worker so far"
                % self.state["spoolDropped"])
        return dropped

    # ---- durable copies ---------------------------------------------------
    def uploadCycle(self, slot):
        """Copy new points, then the checkpoint that follows them, to the store.

        Order matters.  The checkpoint is a hard link snapshot taken first, so
        the points read afterwards include everything flushed before that
        checkpoint was written (the client flushes the dp file, then saves).
        If this process dies between the two uploads the store holds extra
        points and an older checkpoint, which a resume merely re-reports.

        Anything an earlier cycle failed to send goes first, for the same
        reason: points before the checkpoint that follows them."""
        if self.leaseLost:
            raise RuntimeError("lease lost: refusing to publish")
        if not self.slots.heartbeat(slot, self.owner, {}):
            self.leaseLost = True
            self.stopping = True
            raise RuntimeError("lease lost before upload")
        self.lastBeatSuccess = time.monotonic()
        self.drainSpool(slot)
        if self.spoolPending():
            # Cutting a second delta now would spool [offset, more) beside the
            # [offset, less) that just failed -- a superset under a different
            # key, once per cycle, for as long as the outage lasts.  Stopping
            # here keeps the spool one entry deep: those records are still in
            # dp.bin, which rotateDpFile will not remove while the spool holds
            # anything.  It also holds the publication order, since uploading
            # the checkpoint that follows points the store does not have is
            # the one reordering a resume cannot repair.
            raise RuntimeError("store unreachable: %d bytes of points still spooled"
                               % self.spoolBytes())
        snap = self.ckptPath + ".snap"
        if os.path.exists(snap):
            os.remove(snap)
        haveSnap = False
        if os.path.exists(self.ckptPath):
            os.link(self.ckptPath, snap)
            haveSnap = True
        size = os.path.getsize(self.dpPath) if os.path.exists(self.dpPath) else 0
        whole = size - size % RECORD_BYTES
        offset = int(self.state.get("dpOffset", 0))
        if whole > offset:
            delta = os.path.join(self.work, "delta.bin")
            with open(self.dpPath, "rb") as src, open(delta, "wb") as out:
                src.seek(offset)
                out.write(src.read(whole - offset))
            key = "dp/slot-%05d/%s-%016d-%s.bin" % (slot, self.streamId, offset, sha256File(delta))
            metaPath = None
            if self.contract:
                metaPath = delta + ".json"
                writeJson(metaPath, envelope(delta, self.contract, "dp", owner=self.owner, offset=offset))
            entry = self.spoolDelta(delta, key, slot, offset, metaPath)
            self.enforceSpoolBudget()
            try:
                # A raise here leaves the delta spooled and dpOffset where it
                # was; the next cycle drains it before cutting again.
                self.uploadSpoolEntry(entry, slot, source=delta, metaSource=metaPath, check=False)
            finally:
                for path in (delta, delta + ".json"):
                    if os.path.exists(path):
                        os.remove(path)
        if haveSnap:
            it = checkpointIter(snap)
            if it >= 0 and it != self.state.get("ckptIter", -1):
                if self.contract:
                    key = "ckpt/slot-%05d/%s.ck" % (slot, sha256File(snap))
                    self.store.put(snap, key)
                    writeJson(snap + ".json", envelope(snap, self.contract, "checkpoint", iteration=it))
                    self.store.put(snap + ".json", key + ".json")
                    if not self.slots.heartbeat(slot, self.owner, {"checkpointKey": key, "ckptIter": it}):
                        self.leaseLost = True
                        self.stopping = True
                        raise RuntimeError("lease lost: checkpoint pointer not committed")
                    self.lastBeatSuccess = time.monotonic()
                else:
                    self.store.put(snap, self.ckptKey(slot))
                self.state["ckptIter"] = it
                self.saveState()
            os.remove(snap)

    # ---- one client run ---------------------------------------------------
    def runClient(self, slot):
        """Run the client until it exits or a stop/restart is due.

        Returns (returncode, solvedLine)."""
        cmd = self.clientCommand(slot)
        env = dict(os.environ)
        if self.cfg.get("packed", False):
            env["CUDA_VISIBLE_DEVICES"] = str(self.gpu)
        log("starting: " + " ".join(cmd))
        self.proc = subprocess.Popen(cmd, cwd=self.work, env=env, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, text=True, bufsize=1)
        lines = queue.Queue()

        def pump(stream):
            for line in stream:
                lines.put(line.rstrip("\n"))
            lines.put(None)
        threading.Thread(target=pump, args=(self.proc.stdout,), daemon=True).start()

        started = time.time()
        uploadEvery = float(self.cfg.get("uploadEvery", self.cfg["checkpointEvery"]))
        restartAfter = float(self.cfg.get("restartHours", 0)) * 3600
        lastUpload = started
        lastBeat = 0.0
        termSent = None
        restartDue = False
        last = None
        solved = None
        verifiedSolution = False
        tail = []
        eof = False
        while not eof or self.proc.poll() is None:
            try:
                line = lines.get(timeout=1.0)
                if line is None:
                    eof = True
                else:
                    tail = (tail + [line])[-30:]
                    banner = BANNER_RE.search(line)
                    if banner:
                        walks = int(banner.group(1))
                        if walks != int(self.state.get("walks", 0)):
                            self.state["walks"] = walks
                            self.saveState()
                        log("client: " + line)
                    prog = PROGRESS_RE.search(line)
                    if prog:
                        last = {"rate": float(prog.group(2)) * 1e6, "iters": int(prog.group(3)),
                                "dp": int(prog.group(4)), "stored": int(prog.group(5)),
                                "dropped": int(prog.group(6) or 0)}
                    elif not banner:
                        if re.fullmatch(r"\s*k = [0-9]+\s*", line):
                            solved = line.strip()
                            verifiedSolution = False
                        if line.strip() == "verified [k]P == Q" and solved:
                            verifiedSolution = True
                        log("client: " + line)
            except queue.Empty:
                pass
            now = time.time()
            if now - lastBeat >= HEARTBEAT_SECONDS:
                lastBeat = now
                fields = {"ckptIter": int(self.state.get("ckptIter", -1)),
                          "dpUploaded": int(self.state.get("dpUploaded", 0)),
                          "walks": int(self.state.get("walks", 0)),
                          # Unsent and unsendable points, so a worker whose
                          # uploads are failing is visible without its log.
                          "spoolBytes": self.spoolBytes(),
                          "spoolDropped": int(self.state.get("spoolDropped", 0)),
                          "binary": self.cfg.get("binaryKey", ""),
                          "binarySha256": self.cfg.get("binarySha256", ""),
                          "kernelVersion": int(self.cfg.get("kernelVersion") or 0)}
                if last:
                    fields.update(rate=last["rate"], iters=last["iters"], dp=last["dp"], dropped=last["dropped"])
                try:
                    if not self.slots.heartbeat(slot, self.owner, fields):
                        log("lost the lease on slot %d; stopping the client" % slot)
                        self.stopping = True
                        self.leaseLost = True
                    else:
                        self.lastBeatSuccess = time.monotonic()
                except Exception as e:
                    log("heartbeat failed (will retry): %s" % e)
                if self.lastBeatSuccess is None or time.monotonic() - self.lastBeatSuccess >= LEASE_SECONDS:
                    self.leaseLost = True
                    self.stopping = True
                if last:
                    log("%.3f B it/s, %d iterations this run, %d dp, %d uploaded, checkpoint at %d"
                        % (last["rate"] / 1e9, last["iters"], last["dp"],
                           int(self.state.get("dpUploaded", 0)), int(self.state.get("ckptIter", -1))))
                if not restartDue and self.campaignPointerChanged():
                    restartDue = True
                    log("campaign binary moved; checkpointing to pick up the new client")
            if now - lastUpload >= uploadEvery:
                lastUpload = now
                try:
                    self.uploadCycle(slot)
                except Exception as e:
                    log("upload failed (will retry): %s" % e)
            if restartAfter and now - started > restartAfter and not restartDue:
                restartDue = True
                log("scheduled restart after %.1f h" % (restartAfter / 3600))
            if (self.stopping or restartDue) and termSent is None and self.proc.poll() is None:
                termSent = now
                self.proc.send_signal(signal.SIGTERM)
            if termSent and now - termSent > GRACE_SECONDS and self.proc.poll() is None:
                log("client ignored SIGTERM for %d s; killing it" % GRACE_SECONDS)
                self.proc.kill()
        rc = self.proc.wait()
        self.proc = None
        if rc != 0:
            log("client exited %d; last lines:\n  " % rc + "\n  ".join(tail[-12:]))
        else:
            log("client exited 0")
        try:
            self.uploadCycle(slot)
        except Exception as e:
            log("final upload failed: %s" % e)
        return rc, solved if rc == 0 and verifiedSolution and not self.leaseLost else None, restartDue

    def rotateDpFile(self):
        """After the client has exited and every record is uploaded, start a
        fresh dp file so a long-lived instance does not fill its disk.  The
        client reloads its own dp file at startup, which is only the points of
        this stretch, and the campaign's collision detection is the merge's.

        Refused while the spool holds anything: dp.bin is the other local copy
        of those records, and dropping it during an outage would turn a delay
        into a loss."""
        if self.spoolPending():
            log("not rotating dp.bin: %d bytes of points are still unsent" % self.spoolBytes())
            return
        if os.path.exists(self.dpPath):
            size = os.path.getsize(self.dpPath)
            if size - size % RECORD_BYTES <= int(self.state.get("dpOffset", 0)):
                # Reset BEFORE unlink: a crash can re-upload duplicates but
                # cannot skip the prefix of the next file.
                self.state["dpOffset"] = 0
                self.saveState()
                os.remove(self.dpPath)
                self.streamId = uuid.uuid4().hex

    # ---- main loop --------------------------------------------------------
    def run(self):
        self.loadConfig()
        if self.store.exists("solution.json"):
            log("the campaign already has a solution; nothing to do")
            return 0
        failures = 0
        slot = None
        while not self.stopping:
            if slot is None:
                slot = self.claimSlot()
            rc, solved, restartDue = self.runClient(slot)
            if solved:
                res = {"slot": slot, "runId": slot + 1, "owner": self.owner, "line": solved,
                       "when": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
                path = os.path.join(self.work, "solution.json")
                writeJson(path, res)
                self.store.put(path, "solution/slot-%05d.json" % slot)
                self.store.put(path, "solution.json")
                self.slots.release(slot, self.owner, state="solved", extra={"solution": solved})
                log("SOLVED on slot %d: %s" % (slot, solved))
                return 0
            if rc == 6:
                if self.contract:
                    self.slots.release(slot, self.owner, state="error", extra={"reason": "checkpoint rejected; preserve for recovery"})
                    return 1
                self.retireSlot(slot, "checkpoint refused by the client")
                slot = None
                continue
            if rc in (3, 7, 8, 9):
                self.slots.release(slot, self.owner, state="error", extra={"reason": "integrity failure, exit %d" % rc})
                log("reference, report-loss or persistence failure; refusing automatic retries")
                return 1
            if self.stopping:
                break
            if restartDue and rc == 0:
                self.rotateDpFile()
                try:
                    self.reloadClient()
                except Exception as e:
                    log("reload after rollout failed: %s" % e)
                    failures += 1
                    if failures > 5:
                        self.slots.release(slot, self.owner, state="error",
                                           extra={"reason": "rollout reload failed"})
                        return 1
                continue
            failures += 1
            if failures > 5:
                self.slots.release(slot, self.owner, state="error", extra={"reason": "repeated exit %d" % rc})
                log("giving up after repeated failures")
                return 1
            log("client exit %d; retrying in 60 s" % rc)
            for _ in range(60):
                if self.stopping:
                    break
                time.sleep(1)
        if slot is not None:
            try:
                self.slots.release(slot, self.owner)
                log("released slot %d" % slot)
            except Exception as e:
                log("release failed: %s" % e)
        return 0


if __name__ == "__main__":
    try:
        sys.exit(Worker().run())
    except Exception as e:
        log("fatal: %s" % e)
        sys.exit(1)
