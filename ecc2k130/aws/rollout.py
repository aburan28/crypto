#!/usr/bin/env python3
"""Geometry gate and versioned records for a CUDA-kernel rollout.

A compile is not a fleet flip. build.sh --stage publishes bin/<sha>/ and a
rollouts/<sha>.json record. activate assigns the next kernelVersion, writes
an immutable kernels/<n>.json, and rewrites campaign.json only when the
staged prefix keeps the checkpoint shape and the collision contract.

kernelProtocol (ecc2k-kernel-v1) versions the *client pointer*. It is not
storageProtocol and does not migrate the DP corpus. The live unversioned
store can adopt v1 without touching points or checkpoints.

Frozen (activate refuses):
  campaign: curve, dpWeight, workers, batch, blockThreads, minBlocks, walk
  knobs:    BATCH, THREADS, MINBLOCKS, WALK_TABLE,
            PACKED_COMPACT_STATE, PACKED_STATE_TILE
  arches:   staged must be a superset of the live manifest (mixed fleet
            stays on one fat 89+120 binary)
  storage:  a campaign that already has storageProtocol (binary hashes sit
            in that campaign id)

Allowed to move: PACKED_CLMAD and the other ALU / product-pipe knobs,
plus the pointer fields including kernelProtocol / kernelVersion.

max-iters is separate from a rollout: it raises campaign.json maxIters and
touches nothing else (maxItersReasons says when it refuses).

No AWS imports. rollout.sh fetches objects and calls check / apply / adoption
/ max-iters.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time


KERNEL_PROTOCOL = "ecc2k-kernel-v1"
FROZEN_CAMPAIGN = (
    "curve", "dpWeight", "workers", "batch", "blockThreads", "minBlocks", "walk",
)
FROZEN_KNOBS = (
    "BATCH", "THREADS", "MINBLOCKS", "WALK_TABLE",
    "PACKED_COMPACT_STATE", "PACKED_STATE_TILE",
)
WALK_TABLE = {"sigma": "0", "table": "1"}
POINTER_FIELDS = (
    "binaryKey", "hostBinaryKey", "sourceSha256", "binarySha256", "hostBinarySha256",
    "kernelProtocol", "kernelVersion",
)
# First Ada publish on an sm_120-only prefix is bootstrap's job (POINT_CAMPAIGN=1).
# activate of a live mixed fleet must not drop either arch.
MIXED_ARCHES = frozenset({"89", "120"})
SHA256_HEX = re.compile(r"^[0-9a-f]{64}$")


def parseKnobs(text):
    out = {}
    for tok in (text or "").split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def sha256Ok(value):
    return isinstance(value, str) and bool(SHA256_HEX.fullmatch(value))


def nextKernelVersion(campaign):
    return int(campaign.get("kernelVersion") or 0) + 1


def geometryReasons(liveCampaign, stagedManifest, liveManifest=None, allowStrict=False):
    """Human-readable reasons activate must refuse. Empty means the pointer may move."""
    reasons = []
    knobs = parseKnobs(stagedManifest.get("knobs", ""))
    for name in ("BATCH", "THREADS", "MINBLOCKS", "WALK_TABLE"):
        if name not in knobs:
            reasons.append("staged manifest knobs omit %s" % name)
    if reasons:
        return reasons

    campGeom = (
        int(liveCampaign["batch"]),
        int(liveCampaign["blockThreads"]),
        int(liveCampaign["minBlocks"]),
    )
    builtGeom = (int(knobs["BATCH"]), int(knobs["THREADS"]), int(knobs["MINBLOCKS"]))
    if campGeom != builtGeom:
        reasons.append("campaign geometry %r differs from staged knobs %r" % (campGeom, builtGeom))

    walk = liveCampaign.get("walk", "sigma")
    if walk not in WALK_TABLE:
        reasons.append("campaign.json names an unknown walk: %s" % walk)
    elif WALK_TABLE[walk] != knobs["WALK_TABLE"]:
        reasons.append("campaign walk %r differs from staged WALK_TABLE=%s" % (walk, knobs["WALK_TABLE"]))

    liveKnobs = parseKnobs((liveManifest or {}).get("knobs", ""))
    for name in FROZEN_KNOBS:
        if name in liveKnobs and knobs.get(name) != liveKnobs[name]:
            reasons.append("frozen knob %s %s -> %s" % (name, liveKnobs[name], knobs.get(name)))

    stagedArches = {str(a) for a in stagedManifest.get("arches", [])}
    if liveManifest is not None:
        liveArches = {str(a) for a in liveManifest.get("arches", [])}
        missing = sorted(liveArches - stagedArches)
        if missing:
            reasons.append("staged arches %s drop live %s" % (sorted(stagedArches), missing))
    elif liveCampaign.get("binaryKey"):
        missing = sorted(MIXED_ARCHES - stagedArches)
        if missing:
            reasons.append("live campaign has no manifest; staged arches must cover %s (missing %s)"
                           % (sorted(MIXED_ARCHES), missing))

    proto = liveCampaign.get("kernelProtocol") or KERNEL_PROTOCOL
    if proto != KERNEL_PROTOCOL:
        reasons.append("unknown kernelProtocol %r" % proto)
    if liveCampaign.get("storageProtocol") and not allowStrict:
        reasons.append("storageProtocol=%s binds binary hashes into the campaign id; "
                       "kernel rollout needs a new campaign namespace "
                       "(KERNEL_ALLOW_STRICT=1 only after that review)"
                       % liveCampaign["storageProtocol"])
    for key in ("binarySha256", "hostBinarySha256", "sourceSha256"):
        if not sha256Ok(stagedManifest.get(key, "")):
            reasons.append("versioned kernel requires pinned %s" % key)
    return reasons


def applyPointer(campaign, prefix, stagedManifest, kernelVersion=None):
    """Rewrite only the binary pointer fields. Frozen campaign fields stay put."""
    prefix = prefix.rstrip("/")
    out = dict(campaign)
    out["binaryKey"] = prefix + "/ecc2k130"
    out["hostBinaryKey"] = prefix + "/ecc2k130-cpu"
    if stagedManifest.get("sourceSha256"):
        out["sourceSha256"] = stagedManifest["sourceSha256"]
    if stagedManifest.get("binarySha256"):
        out["binarySha256"] = stagedManifest["binarySha256"]
    if stagedManifest.get("hostBinarySha256"):
        out["hostBinarySha256"] = stagedManifest["hostBinarySha256"]
    out["kernelProtocol"] = KERNEL_PROTOCOL
    out["kernelVersion"] = int(kernelVersion) if kernelVersion is not None else nextKernelVersion(campaign)
    return out


def pointerOnly(before, after):
    """True if after differs from before only in POINTER_FIELDS."""
    keys = set(before) | set(after)
    for key in keys:
        if key in POINTER_FIELDS:
            continue
        if before.get(key) != after.get(key):
            return False
    return True


def stagedRecord(prefix, manifest, status="staged"):
    prefix = prefix.rstrip("/")
    return {
        "kernelProtocol": KERNEL_PROTOCOL,
        "prefix": prefix,
        "binaryKey": prefix + "/ecc2k130",
        "hostBinaryKey": prefix + "/ecc2k130-cpu",
        "sourceSha256": manifest.get("sourceSha256", ""),
        "binarySha256": manifest.get("binarySha256", ""),
        "hostBinarySha256": manifest.get("hostBinarySha256", ""),
        "buildSha256": manifest.get("buildSha256", ""),
        "arches": list(manifest.get("arches", [])),
        "knobs": manifest.get("knobs", ""),
        "status": status,
        "stagedAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def kernelRecord(prefix, manifest, version, previous=None, status="active", replaced=""):
    rec = stagedRecord(prefix, manifest, status)
    rec["kernelVersion"] = int(version)
    rec["previousVersion"] = previous
    rec["replaced"] = replaced
    rec["activatedAt"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    return rec


def walkingAdoption(slots, binaryKey, binarySha256="", kernelVersion=None):
    """Count walking leases that have heartbeated the new client.

    A slot is walking when state is active and the lease has not expired.
    Adopted: walking and (binary == binaryKey or binarySha256 matches
    or kernelVersion matches).
    Stale: walking, reported a different binary / version.
    Unknown: walking, no binary or version field (worker.py from before the watch).
    """
    now = int(time.time())
    wantVer = None if kernelVersion in (None, "") else int(kernelVersion)
    walking = adopted = stale = unknown = 0
    for it in slots:
        if it.get("state") not in (None, "active"):
            continue
        if int(it.get("leaseUntil") or 0) < now:
            continue
        walking += 1
        reported = it.get("binary") or ""
        reportedSha = it.get("binarySha256") or ""
        reportedVer = it.get("kernelVersion")
        verHit = wantVer is not None and reportedVer not in (None, "") and int(reportedVer) == wantVer
        if (binaryKey and reported == binaryKey) or (binarySha256 and reportedSha == binarySha256) or verHit:
            adopted += 1
        elif reported or reportedSha or reportedVer not in (None, ""):
            stale += 1
        else:
            unknown += 1
    return {
        "walking": walking,
        "adopted": adopted,
        "stale": stale,
        "unknown": unknown,
        "done": walking > 0 and stale == 0 and unknown == 0,
    }


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


# maxIters restarts a walk that has gone that many steps without a report.
# It is a guard, not part of the walk: a distinguished point is fixed by its
# seed, the walk and dpWeight, and a checkpoint does not record the guard, so
# raising it keeps every point and every checkpoint.  Workers take the new
# value when their client next restarts (restartHours, or a rollout) and
# re-read campaign.json.  Lowering it would cut trails already under way, and
# under storageProtocol the value is hashed into the campaign id, so there a
# raise is a new namespace rather than an edit (WALK-CONSTANT.md section 11.4).
GUARD_PERIOD = 4096           # include/kernel.h ECC_GUARD_PERIOD
MAX_ITERS_LIMIT = (1 << 64) - GUARD_PERIOD   # src/main.cu refuses anything above


def maxItersReasons(campaign, value):
    """Why raising the live campaign's maxIters to `value` is refused, if it is."""
    reasons = []
    current = campaign.get("maxIters")
    if type(value) is not int or not 0 < value <= MAX_ITERS_LIMIT:
        reasons.append("maxIters %r is not a positive step count the client accepts" % (value,))
    if type(current) is not int or current < 0:
        reasons.append("live maxIters %r is unreadable" % (current,))
    elif current == 0:
        reasons.append("live maxIters is 0 (no guard); a positive value would restart walks "
                       "the live store lets run")
    elif type(value) is int and value <= current:
        reasons.append("maxIters %d -> %d is not a raise; only a raise keeps trails already "
                       "under way" % (current, value))
    if campaign.get("storageProtocol"):
        reasons.append("storageProtocol=%s hashes maxIters into the campaign id; a raise there "
                       "is a new campaign namespace, not an edit" % campaign["storageProtocol"])
    return reasons


def setMaxIters(campaign, value):
    """The campaign with maxIters replaced and every other field as it was."""
    out = dict(campaign)
    out["maxIters"] = value
    return out


def _allowStrict(args):
    return bool(getattr(args, "allow_strict", False))


def _cmdCheck(args):
    live = _load(args.campaign)
    staged = _load(args.staged)
    liveMan = _load(args.live_manifest) if args.live_manifest else None
    reasons = geometryReasons(live, staged, liveMan, allowStrict=_allowStrict(args))
    json.dump({"ok": not reasons, "reasons": reasons,
               "nextKernelVersion": nextKernelVersion(live)}, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0 if not reasons else 2


def _cmdApply(args):
    live = _load(args.campaign)
    staged = _load(args.staged)
    liveMan = _load(args.live_manifest) if args.live_manifest else None
    reasons = geometryReasons(live, staged, liveMan, allowStrict=_allowStrict(args))
    if reasons:
        json.dump({"ok": False, "reasons": reasons}, sys.stdout, indent=1)
        sys.stdout.write("\n")
        return 2
    version = int(args.kernel_version) if args.kernel_version else nextKernelVersion(live)
    out = applyPointer(live, args.prefix, staged, kernelVersion=version)
    if not pointerOnly(live, out):
        print("applyPointer touched a frozen field", file=sys.stderr)
        return 1
    with open(args.campaign, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write("\n")
    rec = kernelRecord(args.prefix, staged, version,
                       previous=int(live.get("kernelVersion") or 0) or None,
                       status="active", replaced=live.get("binaryKey") or "")
    if args.record:
        with open(args.record, "w") as fh:
            json.dump(rec, fh, indent=1, sort_keys=True)
            fh.write("\n")
    json.dump({"ok": True, "binaryKey": out["binaryKey"],
               "kernelProtocol": out["kernelProtocol"],
               "kernelVersion": out["kernelVersion"]}, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0


def _cmdAdoption(args):
    slots = _load(args.slots)
    if isinstance(slots, dict):
        slots = [dict(v, slot=int(k)) for k, v in slots.items()]
    out = walkingAdoption(slots, args.binary_key, args.binary_sha256 or "",
                          args.kernel_version)
    json.dump(out, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0 if out["done"] else 3


def _cmdMaxIters(args):
    live = _load(args.campaign)
    reasons = maxItersReasons(live, args.value)
    report = {"ok": not reasons, "reasons": reasons, "from": live.get("maxIters"), "to": args.value}
    if not reasons:
        out = setMaxIters(live, args.value)
        if {k: v for k, v in out.items() if k != "maxIters"} != \
                {k: v for k, v in live.items() if k != "maxIters"}:
            print("setMaxIters touched another field", file=sys.stderr)
            return 1
        with open(args.campaign, "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
            fh.write("\n")
        report["restartHours"] = live.get("restartHours")
    json.dump(report, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0 if not reasons else 2


def _cmdNext(args):
    live = _load(args.campaign)
    json.dump({"kernelProtocol": KERNEL_PROTOCOL,
               "nextKernelVersion": nextKernelVersion(live)}, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    def addGate(sp):
        sp.add_argument("--campaign", required=True)
        sp.add_argument("--staged", required=True, help="staged manifest.json")
        sp.add_argument("--live-manifest")
        sp.add_argument("--allow-strict", action="store_true")

    c = sub.add_parser("check", help="print whether activate may rewrite the pointer")
    addGate(c)
    c.set_defaults(fn=_cmdCheck)

    a = sub.add_parser("apply", help="rewrite campaign.json pointer fields in place")
    addGate(a)
    a.add_argument("--prefix", required=True)
    a.add_argument("--kernel-version", type=int)
    a.add_argument("--record", help="write the immutable kernels/<n>.json body here")
    a.set_defaults(fn=_cmdApply)

    d = sub.add_parser("adoption", help="count walking slots on the new binary")
    d.add_argument("--slots", required=True)
    d.add_argument("--binary-key", required=True)
    d.add_argument("--binary-sha256", default="")
    d.add_argument("--kernel-version", default="")
    d.set_defaults(fn=_cmdAdoption)

    n = sub.add_parser("next-version", help="print the next kernelVersion for a campaign")
    n.add_argument("--campaign", required=True)
    n.set_defaults(fn=_cmdNext)

    m = sub.add_parser("max-iters", help="raise campaign.json maxIters in place, nothing else")
    m.add_argument("--campaign", required=True)
    m.add_argument("--value", required=True, type=int)
    m.set_defaults(fn=_cmdMaxIters)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
