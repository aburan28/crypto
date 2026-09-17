#!/usr/bin/env python3
"""Geometry gate for a CUDA-kernel rollout.

A compile is not a fleet flip. build.sh --stage publishes bin/<sha>/ and a
rollouts/<sha>.json record. activate rewrites campaign.json only when the
staged prefix keeps the checkpoint shape and the collision contract.

Frozen (activate refuses):
  campaign: curve, dpWeight, workers, batch, blockThreads, minBlocks, walk
  knobs:    BATCH, THREADS, MINBLOCKS, WALK_TABLE,
            PACKED_COMPACT_STATE, PACKED_STATE_TILE
  arches:   staged must be a superset of the live manifest (mixed fleet
            stays on one fat 89+120 binary)

Allowed to move: PACKED_CLMAD and the other ALU / product-pipe knobs.

No AWS imports. rollout.sh fetches objects and calls check / apply / adoption.
"""
from __future__ import annotations

import argparse
import json
import sys
import time


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
)
# First Ada publish on an sm_120-only prefix is bootstrap's job (POINT_CAMPAIGN=1).
# activate of a live mixed fleet must not drop either arch.
MIXED_ARCHES = frozenset({"89", "120"})


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


def geometryReasons(liveCampaign, stagedManifest, liveManifest=None):
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

    if not stagedManifest.get("binarySha256"):
        reasons.append("staged manifest has no binarySha256")
    if not liveCampaign.get("binaryKey"):
        reasons.append("campaign.json has no binaryKey; use build.sh (POINT_CAMPAIGN=1) for the first publish")
    return reasons


def applyPointer(campaign, prefix, stagedManifest):
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


def walkingAdoption(slots, binaryKey, binarySha256=""):
    """Count walking leases that have heartbeated the new client.

    A slot is walking when state is active and the lease has not expired.
    Adopted: walking and (binary == binaryKey or binarySha256 matches).
    Stale: walking, reported a different binary.
    Unknown: walking, no binary field (worker.py from before the watch).
    """
    now = int(time.time())
    walking = adopted = stale = unknown = 0
    for it in slots:
        if it.get("state") not in (None, "active"):
            continue
        if int(it.get("leaseUntil") or 0) < now:
            continue
        walking += 1
        reported = it.get("binary") or ""
        reportedSha = it.get("binarySha256") or ""
        if (binaryKey and reported == binaryKey) or (binarySha256 and reportedSha == binarySha256):
            adopted += 1
        elif reported or reportedSha:
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
    return False


def _cmdCheck(args):
    live = _load(args.campaign)
    staged = _load(args.staged)
    liveMan = _load(args.live_manifest) if args.live_manifest else None
    reasons = geometryReasons(live, staged, liveMan)
    json.dump({"ok": not reasons, "reasons": reasons}, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0 if not reasons else 2


def _cmdApply(args):
    live = _load(args.campaign)
    staged = _load(args.staged)
    liveMan = _load(args.live_manifest) if args.live_manifest else None
    reasons = geometryReasons(live, staged, liveMan)
    if reasons:
        json.dump({"ok": False, "reasons": reasons}, sys.stdout, indent=1)
        sys.stdout.write("\n")
        return 2
    out = applyPointer(live, args.prefix, staged)
    if not pointerOnly(live, out):
        print("applyPointer touched a frozen field", file=sys.stderr)
        return 1
    with open(args.campaign, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write("\n")
    json.dump({"ok": True, "binaryKey": out["binaryKey"]}, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0


def _cmdAdoption(args):
    slots = _load(args.slots)
    if isinstance(slots, dict):
        slots = [dict(v, slot=int(k)) for k, v in slots.items()]
    out = walkingAdoption(slots, args.binary_key, args.binary_sha256 or "")
    json.dump(out, sys.stdout, indent=1)
    sys.stdout.write("\n")
    return 0 if out["done"] else 3


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("check", help="print whether activate may rewrite the pointer")
    c.add_argument("--campaign", required=True)
    c.add_argument("--staged", required=True, help="staged manifest.json")
    c.add_argument("--live-manifest")
    c.set_defaults(fn=_cmdCheck)

    a = sub.add_parser("apply", help="rewrite campaign.json pointer fields in place")
    a.add_argument("--campaign", required=True)
    a.add_argument("--staged", required=True)
    a.add_argument("--prefix", required=True)
    a.add_argument("--live-manifest")
    a.set_defaults(fn=_cmdApply)

    d = sub.add_parser("adoption", help="count walking slots on the new binary")
    d.add_argument("--slots", required=True)
    d.add_argument("--binary-key", required=True)
    d.add_argument("--binary-sha256", default="")
    d.set_defaults(fn=_cmdAdoption)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
