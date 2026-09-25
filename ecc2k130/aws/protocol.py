"""Strict seed/orbit corpus protocol; payloads remain 32-byte little-endian records.

Manifests are commit markers, uploaded AFTER immutable payloads. They bind
bytes to the campaign, not to a fictitious (a,b) coefficient representation.
Legacy raw corpora require an explicit legacy merge; never silently adopt one.
"""
import hashlib
import json
import os
import re
import tempfile

PROTOCOL = "ecc2k-seed-orbit-v1"
RECORD_BYTES = 32

# The iteration function is part of the campaign's identity: distinguished
# points from two different walks never collide usefully, so a corpus is
# bound to exactly one.  "sigma" is the equivariant sigma^j + 1 walk every
# corpus so far was collected with; "table" is the additive walk of
# ITERATION-FUNCTION.md §4 (include/tablewalk.h), built with WALK_TABLE=1.
WALKS = {
    "sigma": "sigma^(3+((normal-weight(x)>>1)&7))(R)+R",
    "table": "R+(-1)^eps(R)*sigma^k(R)(T[(normal-weight(x)>>1)&7]);"
             "k=frobenius-phase(x);eps=pivot-coordinate(y);"
             "cycle-rule=advance-h-on-2-4-or-6-step-return-in-Z[tau]-from-last-five-tags",
}


def sha256File(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for data in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def syncDirectory(path):
    fd = os.open(path or ".", os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomicJson(path, obj):
    parent = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix=".commit-", dir=parent)
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=2, sort_keys=True, allow_nan=False)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        syncDirectory(parent)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def campaignContract(config):
    if config.get("storageProtocol") != PROTOCOL:
        raise ValueError("strict storage requires storageProtocol=" + PROTOCOL)
    if config.get("extraArgs"):
        raise ValueError("strict campaigns forbid extraArgs overriding protocol parameters")
    fields = ("curve", "dpWeight", "maxIters", "packed", "workers", "batch",
              "binarySha256", "hostBinarySha256", "sourceSha256")
    c = {k: config[k] for k in fields}
    if type(c["curve"]) is not int or c["curve"] not in (23, 41, 83, 131):
        raise ValueError("strict storage currently supports normal-basis curves 23/41/83/131")
    for key in ("dpWeight", "maxIters", "workers", "batch"):
        if type(c[key]) is not int or c[key] < 0:
            raise ValueError("invalid " + key)
    if not 0 <= c["dpWeight"] <= c["curve"] or not c["workers"] or not c["batch"]:
        raise ValueError("invalid cutoff or worker geometry")
    if type(c["packed"]) is not bool or (c["packed"] and c["curve"] != 131):
        raise ValueError("invalid packed backend")
    for key in ("binarySha256", "hostBinarySha256", "sourceSha256"):
        if not isinstance(c[key], str) or not re.fullmatch("[0-9a-f]{64}", c[key]):
            raise ValueError("strict campaigns require a pinned " + key)
    walk = config.get("walk", "sigma")
    if walk not in WALKS:
        raise ValueError("unknown walk %r; one of %s" % (walk, sorted(WALKS)))
    c.update(protocol=PROTOCOL, recordBytes=32,
             key="min-normal-basis-x-over-frobenius;negation-quotient",
             walk=WALKS[walk],
             seed="run16-walk32-counter16;splitmix64;128-frobenius-terms",
             coefficients="absent;recover-by-seed-replay;verify-kP-equals-Q")
    raw = json.dumps(c, sort_keys=True, separators=(",", ":")).encode()
    return {"id": hashlib.sha256(raw).hexdigest(), "contract": c}


def envelope(path, campaign, kind, **extra):
    size = os.path.getsize(path)
    if kind == "dp" and size % RECORD_BYTES:
        raise ValueError("partial DP record")
    return dict(extra, protocol=PROTOCOL, campaignId=campaign["id"], kind=kind,
                bytes=size, sha256=sha256File(path),
                records=size // RECORD_BYTES if kind == "dp" else None)


def verifyEnvelope(path, manifest, campaign, kind):
    expected = envelope(path, campaign, kind)
    if any(manifest.get(k) != v for k, v in expected.items()):
        raise ValueError("manifest mismatch, corruption or incompatible campaign: " + str(path))


def bindDirectory(root, campaign):
    path = os.path.join(root, "campaign.lock.json")
    if os.path.exists(path):
        with open(path) as fh:
            if json.load(fh) != campaign:
                raise ValueError("directory belongs to a different campaign")
    else:
        # Never bless an old unauthenticated checkpoint/corpus implicitly.
        if any(os.path.exists(os.path.join(root, name)) for name in ("walk.ck", "dp.bin", "state.json")):
            raise ValueError("legacy data requires an explicit audited migration or new directory")
        atomicJson(path, campaign)
