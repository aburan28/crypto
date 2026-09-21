#!/usr/bin/env python3
"""Merge distinguished-point corpora from many workers and solve any collision.

The client writes fixed-width records: the walk seed, then the canonical orbit
representative of the distinguished point (three little-endian 64-bit words).
A v2 corpus leads with an ECC2KDP2 magic and adds the trail length and the
eight branch counts -- the cairn witness -- for 72 bytes a record.
Two records with the same representative and different seeds are a collision,
and the client's own reload path (--load) recomputes both walks and recovers
the logarithm.  This tool only has to find the pair.

The whole ECC2K-130 DP34 corpus is about 2^35.6 records (~1.7 TB), far beyond one
process's hash table, so records are bucketed by a mixed hash of all key words
into bucket files. A pass sorts every bucket that received
new records, drops exact duplicates (a worker resumed from an older checkpoint
re-reports the same points), and reports adjacent equal keys with distinct
seeds.  Buckets are independent, so passes are cheap and incremental.

    merge.py --work WORK --local DIR            ingest DIR/*.bin, detect, solve
    merge.py --work WORK --s3 s3://bucket/dp/   sync from S3 first
    merge.py --work WORK --detect-only          no solve step

Solving runs the host client (ecc2k130-cpu) on a two-record corpus file; it
needs no GPU.  A recovered k is written to WORK/solution.json and, in S3 mode,
copied next to the corpus so every other tool can see the campaign is over.

No type hints, camelCase identifiers (project convention).
"""

import argparse
import fcntl
import json
import os
import re
import subprocess
import sys
import time

import numpy as np
from protocol import (atomicJson, bindDirectory, campaignContract, syncDirectory,
                      verifyEnvelope, sha256File)

RECORD = np.dtype([("seed", "<u8"), ("k0", "<u8"), ("k1", "<u8"), ("k2", "<u8")])
RECORD_BYTES = RECORD.itemsize  # 32

# Corpus v2 carries the cairn witness besides the point: seed, iters, the
# orbit key, and the eight per-branch step counts (see ../CAIRN-WITNESS.md).
# The magic is what tells the two apart -- framing on size alone would read a
# truncated v1 file as v2 and mis-frame every record after the first.
RECORD_V2 = np.dtype([("seed", "<u8"), ("iters", "<u8"),
                      ("k0", "<u8"), ("k1", "<u8"), ("k2", "<u8"),
                      ("counts", "<u4", 8)])
RECORD_V2_BYTES = RECORD_V2.itemsize  # 72
DP_MAGIC_V2 = b"ECC2KDP2"
DP_HEADER_BYTES = 16


def corpusFormat(path):
    """(header bytes, record bytes, dtype) for a corpus file."""
    with open(path, "rb") as fh:
        if fh.read(len(DP_MAGIC_V2)) == DP_MAGIC_V2:
            return DP_HEADER_BYTES, RECORD_V2_BYTES, RECORD_V2
    return 0, RECORD_BYTES, RECORD


def keyRecords(raw, dtype):
    """The (seed, k0, k1, k2) view the merge works on, from either format.

    The witness is deliberately not carried into the buckets.  Merging looks
    for two seeds against one orbit key and nothing else, the solve re-walks
    both trails anyway, and widening every bucket record by 40 bytes to carry
    something the merge never reads would cost the pass its whole margin.
    Whatever wants the witness reads the corpus directly.
    """
    recs = np.frombuffer(raw, dtype=dtype)
    if dtype is RECORD:
        return recs
    out = np.empty(len(recs), dtype=RECORD)
    for field in ("seed", "k0", "k1", "k2"):
        out[field] = recs[field]
    return out


def log(msg):
    # Progress goes to stderr so stdout is exactly one JSON summary.
    print(time.strftime("%H:%M:%S ") + msg, file=sys.stderr, flush=True)


def loadState(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {"version": 2, "buckets": 4096, "offsets": {}, "digests": {}, "bucketHashes": {},
            "dirty": [], "collisions": [], "solved": None}


def saveState(path, state):
    atomicJson(path, state)


def sourceFiles(root):
    out = []
    for dirpath, _, names in os.walk(root):
        for name in names:
            if name.endswith(".bin"):
                out.append(os.path.join(dirpath, name))
    return sorted(out)


def ingest(state, root, work, campaign=None):
    """Append every unread whole record of every source file to its bucket."""
    nb = state["buckets"]
    bucketDir = os.path.join(work, "buckets")
    os.makedirs(bucketDir, exist_ok=True)
    dirty = set(state["dirty"])
    added = 0
    for path in sourceFiles(root):
        rel = os.path.relpath(path, root)
        if campaign:
            metaPath = path + ".json"
            if not os.path.exists(metaPath):
                # A blob without its manifest is not committed yet. Retry next
                # pass; do not advance its offset or count it as accepted work.
                log("pending manifest: " + rel)
                continue
            with open(metaPath) as fh:
                meta = json.load(fh)
            verifyEnvelope(path, meta, campaign, "dp")
            if rel in state["digests"] and state["digests"][rel] != meta["sha256"]:
                raise ValueError("immutable source changed: " + rel)
            state["digests"][rel] = meta["sha256"]
        # Frame by the format the file announces, not by a fixed stride: a v2
        # corpus carries the cairn witness in 72-byte records behind a header,
        # and reading one at 32 would mis-frame every record after the first.
        head, stride, dtype = corpusFormat(path)
        done = max(int(state["offsets"].get(rel, head)), head)
        size = os.path.getsize(path)
        whole = size - (size - head) % stride if size > head else head
        if whole < done:
            raise ValueError("source shrank below committed offset: " + rel)
        if whole == done:
            continue
        with open(path, "rb") as fh:
            fh.seek(done)
            remaining = whole - done
            while remaining:
                data = fh.read(min(remaining, (8 * 1024 * 1024 // stride) * stride))
                if not data or len(data) % stride:
                    raise ValueError("source truncated during ingest: " + rel)
                remaining -= len(data)
                recs = keyRecords(data, dtype)
                # Raw low bits of sparse canonical x keys are NOT uniform.
                h = recs["k0"] * np.uint64(0x9E3779B97F4A7C15)
                h ^= (recs["k1"] + np.uint64(0x632BE59BD9B4E019)) * np.uint64(0xBF58476D1CE4E5B9)
                h ^= (recs["k2"] + np.uint64(0x94D049BB133111EB)) * np.uint64(0x2545F4914F6CDD1D)
                buckets = ((h ^ (h >> np.uint64(29))) & np.uint64(nb - 1)).astype(np.int64)
                order = np.argsort(buckets, kind="stable")
                recs, buckets = recs[order], buckets[order]
                edges = np.flatnonzero(np.diff(buckets)) + 1
                starts, ends = np.concatenate(([0], edges)), np.concatenate((edges, [len(recs)]))
                for s, e in zip(starts, ends):
                    b = int(buckets[s])
                    bucketPath = os.path.join(bucketDir, "%04x.bin" % b)
                    if os.path.exists(bucketPath) and os.path.getsize(bucketPath) % RECORD_BYTES:
                        raise ValueError("torn bucket; rebuild derived merge directory from immutable corpora")
                    with open(bucketPath, "ab") as out:
                        out.write(recs[s:e].tobytes())
                        out.flush()
                        os.fsync(out.fileno())
                    dirty.add(b)
                added += len(recs)
        state["offsets"][rel] = whole
    syncDirectory(bucketDir)
    state["dirty"] = sorted(dirty)
    return added


def detect(state, work):
    """Sort each dirty bucket, drop exact duplicates, return collision pairs."""
    bucketDir = os.path.join(work, "buckets")
    found = []
    total = 0
    for b in list(state["dirty"]):
        path = os.path.join(bucketDir, "%04x.bin" % b)
        if not os.path.exists(path):
            raise ValueError("committed bucket missing: " + path)
        if os.path.getsize(path) % RECORD_BYTES:
            raise ValueError("torn bucket: " + path)
        recs = np.fromfile(path, dtype=RECORD)
        if len(recs) == 0:
            continue
        order = np.lexsort((recs["seed"], recs["k2"], recs["k1"], recs["k0"]))
        recs = recs[order]
        sameKey = (recs["k0"][1:] == recs["k0"][:-1]) & (recs["k1"][1:] == recs["k1"][:-1]) \
            & (recs["k2"][1:] == recs["k2"][:-1])
        sameSeed = recs["seed"][1:] == recs["seed"][:-1]
        keep = np.concatenate(([True], ~(sameKey & sameSeed)))
        recs = recs[keep]
        sameKey = (recs["k0"][1:] == recs["k0"][:-1]) & (recs["k1"][1:] == recs["k1"][:-1]) \
            & (recs["k2"][1:] == recs["k2"][:-1])
        for i in np.flatnonzero(sameKey):
            a, c = recs[i], recs[i + 1]
            found.append({"seedA": "%016x" % int(a["seed"]), "seedB": "%016x" % int(c["seed"]),
                          "key": "%016x%016x%016x" % (int(a["k2"]), int(a["k1"]), int(a["k0"])),
                          "bucket": b})
        tmp = path + ".tmp"
        with open(tmp, "wb") as fh:
            recs.tofile(fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        syncDirectory(bucketDir)
        total += len(recs)
    state["dirty"] = []
    return found, total


def bucketTotal(work):
    bucketDir = os.path.join(work, "buckets")
    if not os.path.isdir(bucketDir):
        return 0
    return sum(os.path.getsize(os.path.join(bucketDir, f)) for f in os.listdir(bucketDir)
               if f.endswith(".bin")) // RECORD_BYTES


def verifyBuckets(state, work):
    for name, digest in state.get("bucketHashes", {}).items():
        path = os.path.join(work, "buckets", name)
        if not os.path.isfile(path) or sha256File(path) != digest:
            raise ValueError("bucket integrity failure; rebuild merge state from immutable source corpora: " + name)


def hashBuckets(state, work):
    directory = os.path.join(work, "buckets")
    state["bucketHashes"] = {name: sha256File(os.path.join(directory, name))
                             for name in os.listdir(directory) if name.endswith(".bin")}


def solve(pair, client, curve, work, extra, timeout=3600):
    """Hand the two colliding records to the host client and parse its answer."""
    seedA, seedB = int(pair["seedA"], 16), int(pair["seedB"], 16)
    key = int(pair["key"], 16)
    words = [(key >> (64 * i)) & ((1 << 64) - 1) for i in range(3)]
    recs = np.array([(seedA,) + tuple(words), (seedB,) + tuple(words)], dtype=RECORD)
    pairFile = os.path.join(work, "pair-%s-%s.bin" % (pair["seedA"], pair["seedB"]))
    recs.tofile(pairFile)
    cmd = [client, "--curve", str(curve), "--threads", "1", "--steps", "1", "--launches", "1",
           "--verify", "0", "--run-id", "65535", "--load", pairFile] + extra
    log("solving: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    out = proc.stdout + proc.stderr
    result = {"pair": pair, "returncode": proc.returncode, "k": None, "verified": False,
              "matchesPublished": None, "tail": out.strip().splitlines()[-12:]}
    for line in out.splitlines():
        t = line.strip()
        if t.startswith("k = "):
            result["k"] = t[4:].strip()
        if "verified [k]P == Q" in t:
            result["verified"] = True
        m = re.search(r"matches the (?:published|planted) (?:solution|discrete log): (\w+)", t)
        if m:
            result["matchesPublished"] = m.group(1)
    result["verified"] = bool(result["verified"] and result["k"] and
                              proc.returncode == 0 and result["matchesPublished"] != "NO")
    return result


def s3Sync(uri, dest):
    os.makedirs(dest, exist_ok=True)
    cmd = ["aws", "s3", "sync", uri, dest, "--only-show-errors"]
    log(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, help="state directory (buckets, offsets, results)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--local", help="directory of *.bin corpus files")
    src.add_argument("--s3", help="s3://bucket/prefix/ of corpus deltas to sync")
    ap.add_argument("--client", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "..", "ecc2k130-cpu"),
                    help="host client used to rewalk and solve (default ../ecc2k130-cpu)")
    ap.add_argument("--curve", type=int, default=131)
    ap.add_argument("--buckets", type=int, default=4096, help="bucket count for a new work dir")
    ap.add_argument("--detect-only", action="store_true")
    ap.add_argument("--campaign", help="strict campaign.json; require committed checksummed chunks")
    ap.add_argument("--legacy", action="store_true", help="explicitly allow unversioned raw corpora (not certified)")
    ap.add_argument("--solve-timeout", type=float, default=3600)
    ap.add_argument("--dp-weight", type=int, default=None,
                    help="the campaign's cutoff; the rewalk must stop at the same points the workers reported")
    ap.add_argument("--client-arg", action="append", default=[],
                    help="extra client argument for the solve step (repeatable)")
    args = ap.parse_args()
    if bool(args.campaign) == bool(args.legacy):
        ap.error("choose exactly one of --campaign or --legacy")
    if args.buckets < 1 or args.buckets & (args.buckets - 1):
        ap.error("--buckets must be a positive power of two")
    if args.solve_timeout <= 0:
        ap.error("--solve-timeout must be positive")
    campaign = None
    if args.campaign:
        if args.client_arg:
            ap.error("strict merge forbids --client-arg overrides")
        with open(args.campaign) as fh:
            config = json.load(fh)
        campaign = campaignContract(config)
        if sha256File(args.client) != config["hostBinarySha256"]:
            ap.error("solver binary hash differs from pinned campaign hostBinarySha256")
        args.curve, args.dp_weight = config["curve"], config["dpWeight"]
        if config["maxIters"]:
            args.client_arg = ["--max-iters", str(config["maxIters"])]
    if args.dp_weight is not None:
        args.client_arg = ["--dp-weight", str(args.dp_weight)] + args.client_arg

    os.makedirs(args.work, exist_ok=True)
    with open(os.path.join(args.work, "merge.lock"), "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return runMerge(args, campaign)


def runMerge(args, campaign):
    if campaign:
        bindDirectory(args.work, campaign)
    statePath = os.path.join(args.work, "state.json")
    state = loadState(statePath)
    if state.get("version") != 2:
        raise ValueError("old bucket layout; use a fresh merge directory (preserve source corpora)")
    verifyBuckets(state, args.work)
    binding = {"campaignId": campaign["id"] if campaign else None,
               "curve": args.curve, "dpWeight": args.dp_weight, "clientArgs": args.client_arg}
    if "binding" in state and state["binding"] != binding:
        raise ValueError("merge state belongs to different campaign/solve parameters")
    state["binding"] = binding
    if not os.path.exists(statePath):
        state["buckets"] = args.buckets
    if state.get("solved"):
        log("already solved: k = %s" % state["solved"]["k"])
        print(json.dumps(state["solved"], indent=1))
        return 0

    root = args.local
    if args.s3:
        root = os.path.join(args.work, "s3cache")
        s3Sync(args.s3, root)
    t0 = time.time()
    added = ingest(state, root, args.work, campaign)
    hashBuckets(state, args.work)
    saveState(statePath, state)
    log("ingested %d new records in %.1f s (%d buckets to re-sort)"
        % (added, time.time() - t0, len(state["dirty"])))
    t1 = time.time()
    found, sorted_ = detect(state, args.work)
    known = {(c["seedA"], c["seedB"]) for c in state["collisions"]}
    new = [c for c in found if (c["seedA"], c["seedB"]) not in known]
    state["collisions"].extend(new)
    hashBuckets(state, args.work)
    saveState(statePath, state)
    # Every recorded pair without a solve result is still pending: a
    # detect-only pass records pairs, a later pass solves them.
    pending = [c for c in state["collisions"] if not (c.get("result") or {}).get("verified")]
    log("sorted %d records in %.1f s; corpus %d records; %d collision(s), %d new, %d unsolved"
        % (sorted_, time.time() - t1, bucketTotal(args.work), len(state["collisions"]), len(new), len(pending)))
    summary = {"added": added, "corpus": bucketTotal(args.work),
               "collisions": len(state["collisions"]), "new": new, "pending": len(pending), "solution": None}
    if args.detect_only or not pending:
        print(json.dumps(summary, indent=1))
        return 0

    for pair in pending:
        try:
            res = solve(pair, args.client, args.curve, args.work, args.client_arg, args.solve_timeout)
        except (OSError, subprocess.TimeoutExpired) as exc:
            pair["result"] = {"verified": False, "error": type(exc).__name__}
            saveState(statePath, state)
            continue
        pair["result"] = {"k": res["k"], "verified": res["verified"], "returncode": res["returncode"],
                          "tail": res["tail"]}
        saveState(statePath, state)
        log("seeds %s / %s: k = %s, verified %s" % (pair["seedA"], pair["seedB"], res["k"], res["verified"]))
        if res["k"] and res["verified"]:
            res["when"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            state["solved"] = res
            saveState(statePath, state)
            saveState(os.path.join(args.work, "solution.json"), res)
            if args.s3:
                dest = args.s3.rstrip("/").rsplit("/", 1)[0] + "/solution.json"
                subprocess.run(["aws", "s3", "cp", os.path.join(args.work, "solution.json"), dest,
                                "--only-show-errors"])
            summary["solution"] = res
            print(json.dumps(summary, indent=1))
            return 0
    print(json.dumps(summary, indent=1))
    return 2


if __name__ == "__main__":
    sys.exit(main())
