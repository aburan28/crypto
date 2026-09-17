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

PROGRESS_RE = re.compile(
    r"([\d.]+)\s+s\s+([\d.]+)\s+M it/s\s+(\d+)\s+iterations\s+(\d+)\s+dp\s+(\d+)\s+stored"
    r"(?:\s+(\d+)\s+dropped)?")
RECORD_BYTES = 32
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


def isCpuDevice():
    return os.environ.get("ECC_DEVICE", "").lower() == "cpu"


def isCpuSlotRecord(item):
    name = str(item.get("gpuName") or "")
    return name.startswith("cpu") or name == "cpu" or "/cpu" in name


def idleSlotClaimable(item, info, newOnly):
    """Whether this claimant may resume an expired/idle slot.

    CPU and GPU share one registry. Mixing a host-shaped checkpoint with the
    packed client (or the reverse) is exit 6, and on the unversioned store
    retireSlot then drops that run id. CPU checkpoints also require an exact
    thread match: cpu/4 must not resume cpu/2. An optional @host suffix on
    gpuName is ignored when comparing.
    """
    cpuSlot = isCpuSlotRecord(item)
    if newOnly:
        have = str(item.get("gpuName") or "").split("@", 1)[0]
        want = str(info.get("gpuName") or "").split("@", 1)[0]
        return cpuSlot and bool(want) and have == want
    return not cpuSlot


def claimNewOnly():
    # CPU (and FPGA host) walks must never resume a GPU checkpoint: loading a
    # packed-shape walk.ck into ecc2k130-cpu exits 6 and the supervisor would
    # retire the slot. Prefer an explicit ECC_CLAIM_NEW=1; CPU mode implies it.
    if os.environ.get("ECC_CLAIM_NEW", "").strip() in ("1", "true", "yes"):
        return True
    return isCpuDevice()


def cpuThreadCount():
    raw = os.environ.get("ECC_THREADS", "").strip()
    if raw:
        return max(1, int(raw))
    return max(1, int(os.cpu_count() or 1))


def gpuName(gpu):
    # Non-GPU clients (CPU / FPGA host) have no nvidia-smi; the bootstrap
    # tells us what the device is instead.
    if os.environ.get("ECC_DEVICE_NAME"):
        return os.environ["ECC_DEVICE_NAME"]
    if isCpuDevice():
        return "cpu/%d" % cpuThreadCount()
    try:
        r = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader", "-i", str(gpu)],
                           capture_output=True, text=True)
    except OSError:
        return "cpu"
    return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else "cpu"


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
                      "--projection-expression", "#s, leaseUntil, #st, #o, gpuName",
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

    def claim(self, owner, info, newOnly=False):
        now = int(time.time())
        items = self.scan()
        free = sorted(it["slot"] for it in items
                      if it.get("state") not in ("retired", "solved", "error")
                      and int(it.get("leaseUntil") or 0) < now
                      and idleSlotClaimable(it, info, newOnly))
        names = {"#o": "owner", "#st": "state"}
        for slot in free:
            values = {":me": dv(owner), ":t": dv(now + LEASE_SECONDS), ":now": dv(now), ":active": dv("active"),
                      ":idle": dv("idle"), ":inst": dv(info["instance"]), ":gpu": dv(info["gpu"]),
                      ":gpuName": dv(info["gpuName"])}
            # Only an expired or released lease may be taken, and only from a
            # slot that is still walking: retired, solved and error slots keep
            # their run id forever so its seeds are never walked twice.
            ok = self._update(slot, "SET #o = :me, leaseUntil = :t, claimedAt = :now, #st = :active, "
                              "instance = :inst, gpu = :gpu, gpuName = :gpuName", names, values,
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
                    "instance": dv(info["instance"]), "gpu": dv(info["gpu"]), "gpuName": dv(info["gpuName"])}
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

    def claim(self, owner, info, newOnly=False):
        now = int(time.time())
        items = self.scan()
        for it in sorted(items, key=lambda x: x["slot"]):
            if it.get("state") not in (None, "active", "idle") or int(it.get("leaseUntil") or 0) >= now:
                continue
            if not idleSlotClaimable(it, info, newOnly):
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

    def claim(self, owner, info, newOnly=False):
        def fn(items):
            now = int(time.time())
            for k in sorted(items, key=int):
                it = items[k]
                if it.get("state") in ("retired", "solved", "error") or int(it.get("leaseUntil") or 0) >= now:
                    continue
                if not idleSlotClaimable(it, info, newOnly):
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
        self.cpu = isCpuDevice()
        defaultClient = "ecc2k130-cpu" if self.cpu else "ecc2k130"
        self.client = os.environ.get("ECC_CLIENT", os.path.join(self.root, defaultClient))
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
        deviceTag = "cpu%d" % self.gpu if self.cpu else "gpu%d" % self.gpu
        self.owner = "%s:%s:%s" % (self.instance, deviceTag, uuid.uuid4().hex)
        self.work = os.path.join(self.root, deviceTag)
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

    def expectedClientSha(self):
        if self.cpu:
            return self.cfg.get("hostBinarySha256") or self.cfg.get("binarySha256")
        return self.cfg.get("binarySha256")

    # ---- campaign configuration ------------------------------------------
    def loadConfig(self):
        path = os.path.join(self.work, "campaign.json")
        if not self.store.get("campaign.json", path):
            raise RuntimeError("campaign.json missing from the store")
        self.cfg = readJson(path)
        for key in ("curve", "steps", "checkpointEvery"):
            if key not in self.cfg:
                raise RuntimeError("campaign.json lacks %r" % key)
        expected = self.expectedClientSha()
        if expected and sha256File(self.client) != expected:
            raise RuntimeError("client binary hash differs from campaign")
        if self.cfg.get("storageProtocol"):
            self.contract = campaignContract(self.cfg)
            bindDirectory(self.work, self.contract)
        elif not os.environ.get("ECC_ALLOW_LEGACY_STORAGE"):
            raise RuntimeError("unversioned campaign: set storageProtocol; legacy storage requires ECC_ALLOW_LEGACY_STORAGE=1")

    def clientCommand(self, slot):
        c = self.cfg
        cmd = [self.client, "--curve", str(c["curve"]), "--steps", str(c["steps"]), "--launches", "0",
               "--run-id", str(slot + 1), "--dp-file", self.dpPath, "--checkpoint", self.ckptPath,
               "--checkpoint-every", str(int(c["checkpointEvery"])), "--verify", str(int(c.get("verify", 0)))]
        if self.cpu:
            # Host binary: never --packed. Thread count is the box, not the
            # GPU preset (campaign workers=385024 would OOM a c7i).
            cmd += ["--threads", str(cpuThreadCount())]
        else:
            if c.get("packed", False):
                cmd += ["--packed", "--device", "0"]
            if c.get("workers"):
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
        info = {"instance": self.instance, "gpu": self.gpu, "gpuName": self.gpuName}
        newOnly = claimNewOnly()
        if newOnly:
            log("CPU-safe claim: resume idle %s slots only; otherwise allocate a new slot" % self.gpuName)
        slot = self.slots.claim(self.owner, info, newOnly=newOnly)
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
            # A different slot than this work dir last held: nothing local applies.
            for name in ("dp.bin", "walk.ck", "walk.ck.snap", "walk.ck.remote"):
                p = os.path.join(self.work, name)
                if os.path.exists(p):
                    os.remove(p)
            self.state = {"slot": slot, "dpOffset": 0, "ckptIter": -1, "dpUploaded": 0}
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
        return slot

    def retireSlot(self, slot, reason):
        log("retiring slot %d: %s" % (slot, reason))
        if os.path.exists(self.ckptPath):
            self.store.put(self.ckptPath, "ckpt/retired/slot-%05d.ck" % slot)
        self.slots.release(slot, self.owner, state="retired", extra={"reason": reason})
        self.state = {}
        self.saveState()

    # ---- durable copies ---------------------------------------------------
    def uploadCycle(self, slot):
        """Copy new points, then the checkpoint that follows them, to the store.

        Order matters.  The checkpoint is a hard link snapshot taken first, so
        the points read afterwards include everything flushed before that
        checkpoint was written (the client flushes the dp file, then saves).
        If this process dies between the two uploads the store holds extra
        points and an older checkpoint, which a resume merely re-reports."""
        if self.leaseLost:
            raise RuntimeError("lease lost: refusing to publish")
        if not self.slots.heartbeat(slot, self.owner, {}):
            self.leaseLost = True
            self.stopping = True
            raise RuntimeError("lease lost before upload")
        self.lastBeatSuccess = time.monotonic()
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
            self.store.put(delta, key)
            if self.contract:
                meta = envelope(delta, self.contract, "dp", owner=self.owner, offset=offset)
                writeJson(delta + ".json", meta)
                self.store.put(delta + ".json", key + ".json")
            os.remove(delta)
            self.state["dpOffset"] = whole
            # Cumulative across dp file rotations, so the dashboard's count
            # is this slot's whole contribution.
            self.state["dpUploaded"] = int(self.state.get("dpUploaded", 0)) + (whole - offset) // RECORD_BYTES
            self.saveState()
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
        if (not self.cpu) and self.cfg.get("packed", False):
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
                    prog = PROGRESS_RE.search(line)
                    if prog:
                        last = {"rate": float(prog.group(2)) * 1e6, "iters": int(prog.group(3)),
                                "dp": int(prog.group(4)), "stored": int(prog.group(5)),
                                "dropped": int(prog.group(6) or 0)}
                    else:
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
                          "binary": self.cfg.get("binaryKey", "")}
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
        this stretch, and the campaign's collision detection is the merge's."""
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
