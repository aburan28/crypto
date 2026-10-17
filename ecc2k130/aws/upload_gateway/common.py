"""Shared validation for the ECC2K-130 upload broker and validator.

This module deliberately uses only the Python standard library.  It is loaded
inside both Lambda functions and by the unit suite without an AWS SDK shim.
"""

import hashlib
import json
import os
import re
import struct
from functools import lru_cache


DP_V2_MAGIC = b"ECC2KDP2"
DP_TABLE3_MAGIC = b"ECC2KDT3"
DP_HEADER_BYTES = 16
DP_V1_BYTES = 32
DP_V2_BYTES = 72
MAX_DP_BYTES = 256 * 1024 * 1024
MAX_CHECKPOINT_BYTES = 768 * 1024 * 1024
MAX_CONTROL_BYTES = 1024 * 1024

DP_KEY_RE = re.compile(
    r"^dp/slot-(\d{5})/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$")
CHECKPOINT_KEY_RE = re.compile(
    r"^ckpt/slot-(\d{5})/([0-9a-f]{64})\.ck$")
LEGACY_CHECKPOINT_KEY_RE = re.compile(r"^ckpt/slot-(\d{5})\.ck$")
RETIRED_CHECKPOINT_KEY_RE = re.compile(r"^ckpt/retired/slot-(\d{5})\.ck$")
SOLUTION_KEY_RE = re.compile(r"^solution/slot-(\d{5})\.json$")
SLOT_KEY_RE = re.compile(r"^slots/slot-(\d{5})\.json$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
STREAM_RE = re.compile(r"^[0-9a-f]{32}$")
WALKS = {
    "sigma": "sigma^(3+((normal-weight(x)>>1)&7))(R)+R",
    "table": "R+(-1)^eps(R)*sigma^k(R)(T[(normal-weight(x)>>1)&7]);"
             "k=frobenius-phase(x);eps=pivot-coordinate(y);"
             "cycle-rule=v3-raw-cycle8-canonical-eligible-anchor",
}


class ValidationError(ValueError):
    pass


def canonicalJson(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode()


def campaignIdentity(config):
    """Stable identity for strict and pre-protocol campaigns.

    Strict campaigns already publish the contract id in their manifests.  The
    live campaign predates that field, so bind its walk-defining inputs here;
    mutable rollout and restart settings intentionally do not split a corpus.
    """
    if config.get("storageProtocol") == "ecc2k-seed-orbit-v1":
        fields = ("curve", "dpWeight", "maxIters", "packed", "workers", "batch",
                  "binarySha256", "hostBinarySha256", "sourceSha256")
        body = {key: config[key] for key in fields}
        body["maxIters"] = config.get("contractMaxIters", body["maxIters"])
        walk = config.get("walk", "sigma")
        if walk not in WALKS:
            raise ValidationError("campaign names an unknown walk")
        body.update(protocol="ecc2k-seed-orbit-v1", recordBytes=32,
                    key="min-normal-basis-x-over-frobenius;negation-quotient",
                    walk=WALKS[walk],
                    seed="run16-walk32-counter16;splitmix64;128-frobenius-terms",
                    coefficients="absent;recover-by-seed-replay;verify-kP-equals-Q")
        return hashlib.sha256(canonicalJson(body)).hexdigest()
    fields = ("curve", "dpWeight", "packed", "workers", "batch", "walk",
              "sourceSha256")
    body = {key: config.get(key) for key in fields}
    body["walk"] = body.get("walk") or "sigma"
    body["storageProtocol"] = config.get("storageProtocol") or "legacy"
    return hashlib.sha256(canonicalJson(body)).hexdigest()


def objectSpec(key):
    match = DP_KEY_RE.fullmatch(key)
    if match:
        return {"kind": "dp", "slot": int(match.group(1)),
                "streamId": match.group(2), "offset": int(match.group(3)),
                "keySha256": match.group(4), "maxBytes": MAX_DP_BYTES}
    match = CHECKPOINT_KEY_RE.fullmatch(key)
    if match:
        return {"kind": "checkpoint", "slot": int(match.group(1)),
                "keySha256": match.group(2), "maxBytes": MAX_CHECKPOINT_BYTES}
    match = LEGACY_CHECKPOINT_KEY_RE.fullmatch(key)
    if match:
        return {"kind": "checkpoint", "slot": int(match.group(1)),
                "maxBytes": MAX_CHECKPOINT_BYTES}
    match = RETIRED_CHECKPOINT_KEY_RE.fullmatch(key)
    if match:
        return {"kind": "retired-checkpoint", "slot": int(match.group(1)),
                "maxBytes": MAX_CHECKPOINT_BYTES}
    match = SOLUTION_KEY_RE.fullmatch(key)
    if match:
        return {"kind": "solution", "slot": int(match.group(1)),
                "maxBytes": MAX_CONTROL_BYTES}
    raise ValidationError("object key is outside the worker upload allowlist")


def validateTicket(body):
    required = ("key", "bytes", "sha256", "workerId", "owner")
    missing = [key for key in required if key not in body]
    if missing:
        raise ValidationError("ticket lacks " + ", ".join(missing))
    spec = objectSpec(body["key"])
    if isinstance(body["bytes"], bool) or not isinstance(body["bytes"], int):
        raise ValidationError("bytes must be an integer")
    if body["bytes"] <= 0 or body["bytes"] > spec["maxBytes"]:
        raise ValidationError("object size is outside the allowed range")
    if not isinstance(body["sha256"], str) or not SHA256_RE.fullmatch(body["sha256"]):
        raise ValidationError("sha256 must be lowercase hexadecimal")
    if spec.get("keySha256") and body["sha256"] != spec["keySha256"]:
        raise ValidationError("object hash does not match its immutable key")
    for field, limit in (("workerId", 128), ("owner", 256)):
        if not isinstance(body[field], str) or not body[field] or len(body[field]) > limit:
            raise ValidationError(field + " is invalid")
    manifest = body.get("manifest")
    if manifest is not None:
        if not isinstance(manifest, dict):
            raise ValidationError("manifest must be an object")
        if manifest.get("bytes") != body["bytes"] or manifest.get("sha256") != body["sha256"]:
            raise ValidationError("manifest does not bind the proposed body")
        if manifest.get("kind") not in (None, spec["kind"]):
            raise ValidationError("manifest kind does not match the object key")
        if manifest.get("owner") not in (None, body["owner"]):
            raise ValidationError("manifest owner does not match the lease owner")
        if spec["kind"] == "dp" and manifest.get("offset") not in (None, spec["offset"]):
            raise ValidationError("manifest offset does not match the object key")
    return spec


@lru_cache(maxsize=None)
def _frobeniusMap(bits):
    """Coordinate permutation used by the client's type-II ONB.

    The stored limbs are the client's *permuted* normal-basis coordinates,
    not the plain cyclic coordinates used in the protocol description.  In
    ScalarOnb::sigma, coordinate i maps to fold(2i) modulo 2m+1.  Applying a
    simple bit rotation therefore rejects valid production records.
    """
    ring = 2 * bits + 1
    result = []
    for source in range(1, bits + 1):
        target = (2 * source) % ring
        if target > bits:
            target = ring - target
        result.append(target - 1)
    return tuple(result)


def _frobenius(value, bits):
    mapping = _frobeniusMap(bits)
    result = 0
    while value:
        lowest = value & -value
        source = lowest.bit_length() - 1
        result |= 1 << mapping[source]
        value ^= lowest
    return result


def isCanonicalOrbit(value, bits):
    current = value
    best = value
    for _ in range(1, bits):
        current = _frobenius(current, bits)
        if current < best:
            best = current
    return value == best


def dpLayout(header):
    if header[:8] == DP_V2_MAGIC:
        if len(header) < DP_HEADER_BYTES:
            raise ValidationError("truncated v2 DP header")
        version, stride = struct.unpack("<II", header[8:16])
        if version != 2 or stride != DP_V2_BYTES:
            raise ValidationError("invalid v2 DP header")
        return DP_HEADER_BYTES, DP_V2_BYTES, 16, True
    if header[:8] == DP_TABLE3_MAGIC:
        if len(header) < DP_HEADER_BYTES:
            raise ValidationError("truncated table-v3 DP header")
        version, stride = struct.unpack("<II", header[8:16])
        if version != 3 or stride != DP_V1_BYTES:
            raise ValidationError("invalid table-v3 DP header")
        return DP_HEADER_BYTES, DP_V1_BYTES, 8, False
    return 0, DP_V1_BYTES, 8, False


def validateDpFile(path, slot, config):
    curve = config.get("curve")
    cutoff = config.get("dpWeight")
    workers = config.get("workers")
    batch = config.get("batch")
    if any(isinstance(value, bool) or not isinstance(value, int)
           for value in (curve, cutoff, workers, batch)):
        raise ValidationError("campaign walk geometry is invalid")
    if curve <= 128 or curve > 192 or not 0 <= cutoff <= curve or workers <= 0 or batch <= 0:
        raise ValidationError("campaign walk geometry is outside the supported range")
    size = os.path.getsize(path)
    with open(path, "rb") as source:
        header = source.read(DP_HEADER_BYTES)
        base, stride, canonOffset, witness = dpLayout(header)
        if size <= base or (size - base) % stride:
            raise ValidationError("DP body is not a whole number of records")
        source.seek(base)
        records = 0
        seeds = set()
        runId = slot + 1
        laneLimit = workers * batch
        topBits = curve - 128
        while True:
            record = source.read(stride)
            if not record:
                break
            if len(record) != stride:
                raise ValidationError("truncated DP record")
            seed = struct.unpack_from("<Q", record)[0]
            if seed >> 48 != runId:
                raise ValidationError("DP seed is outside the slot run id")
            if ((seed >> 16) & 0xffffffff) >= laneLimit:
                raise ValidationError("DP seed is outside the campaign lane range")
            if seed in seeds:
                raise ValidationError("DP body repeats a seed")
            seeds.add(seed)
            k0, k1, k2 = struct.unpack_from("<QQQ", record, canonOffset)
            if k2 >> topBits:
                raise ValidationError("DP representative has bits outside the field")
            point = k0 | (k1 << 64) | (k2 << 128)
            if point.bit_count() > cutoff:
                raise ValidationError("record does not satisfy the campaign DP predicate")
            if not isCanonicalOrbit(point, curve):
                raise ValidationError("record is not the canonical Frobenius representative")
            if witness:
                iterations = struct.unpack_from("<Q", record, 8)[0]
                counts = struct.unpack_from("<8I", record, 40)
                if sum(counts) != iterations:
                    raise ValidationError("v2 witness counts do not sum to its iteration count")
            records += 1
    return {"records": records, "recordBytes": (size - base),
            "headerBytes": base, "stride": stride}


def reserveStreamRange(start, recordBytes, headerBytes, nextOffset=None):
    """Validate and return a DP stream's next source-file byte offset."""
    expected = int(headerBytes) if nextOffset is None else int(nextOffset)
    if int(start) != expected:
        if nextOffset is None:
            raise ValidationError("a new DP stream must begin at its first record")
        raise ValidationError("DP stream offset overlaps or leaves a gap")
    return int(start) + int(recordBytes)


def validateManifest(manifest, ticket, campaignId, records):
    if manifest is None:
        return
    for field, expected in (("bytes", ticket["bytes"]),
                            ("sha256", ticket["sha256"]),
                            ("campaignId", campaignId)):
        if field in manifest and manifest[field] != expected:
            raise ValidationError("manifest %s does not match the validated object" % field)
    if ticket["kind"] == "dp" and "records" in manifest:
        # Historical v1 manifests count 32-byte records.  New broker metadata
        # reports the actual framed record count, so accept either exact form.
        legacy = ticket["bytes"] // DP_V1_BYTES
        if manifest["records"] not in (records, legacy):
            raise ValidationError("manifest record count does not match the body")
