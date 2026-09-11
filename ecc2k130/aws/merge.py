#!/usr/bin/env python3
"""Merge distinguished-point corpora from many workers and solve any collision.

The client writes 32-byte records: the walk seed, then the canonical orbit
representative of the distinguished point (three little-endian 64-bit words).
Two records with the same representative and different seeds are a collision,
and the client's own reload path (--load) recomputes both walks and recovers
the logarithm.  This tool only has to find the pair.

The whole ECC2K-130 corpus is about 2^35.6 records (~425 GB), far beyond one
process's hash table, so records are bucketed by the low bits of the first
key word into fixed-size bucket files.  A pass sorts every bucket that received
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
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

RECORD = np.dtype([("seed", "<u8"), ("k0", "<u8"), ("k1", "<u8"), ("k2", "<u8")])
RECORD_BYTES = RECORD.itemsize  # 32


def log(msg):
    # Progress goes to stderr so stdout is exactly one JSON summary.
    print(time.strftime("%H:%M:%S ") + msg, file=sys.stderr, flush=True)


def loadState(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {"buckets": 4096, "offsets": {}, "dirty": [], "collisions": [], "solved": None}


def saveState(path, state):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=1, sort_keys=True)
    os.replace(tmp, path)


def sourceFiles(root):
    out = []
    for dirpath, _, names in os.walk(root):
        for name in names:
            if name.endswith(".bin"):
                out.append(os.path.join(dirpath, name))
    return sorted(out)


def ingest(state, root, work):
    """Append every unread whole record of every source file to its bucket."""
    nb = state["buckets"]
    bucketDir = os.path.join(work, "buckets")
    os.makedirs(bucketDir, exist_ok=True)
    dirty = set(state["dirty"])
    added = 0
    for path in sourceFiles(root):
        rel = os.path.relpath(path, root)
        done = int(state["offsets"].get(rel, 0))
        size = os.path.getsize(path)
        whole = size - size % RECORD_BYTES
        if whole <= done:
            continue
        with open(path, "rb") as fh:
            fh.seek(done)
            recs = np.frombuffer(fh.read(whole - done), dtype=RECORD)
        buckets = (recs["k0"] & (nb - 1)).astype(np.int64)
        order = np.argsort(buckets, kind="stable")
        recs = recs[order]
        buckets = buckets[order]
        edges = np.flatnonzero(np.diff(buckets)) + 1
        starts = np.concatenate(([0], edges))
        ends = np.concatenate((edges, [len(recs)]))
        for s, e in zip(starts, ends):
            b = int(buckets[s])
            with open(os.path.join(bucketDir, "%04x.bin" % b), "ab") as out:
                out.write(recs[s:e].tobytes())
            dirty.add(b)
        state["offsets"][rel] = whole
        added += len(recs)
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
            continue
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
        recs.tofile(tmp)
        os.replace(tmp, path)
        total += len(recs)
    state["dirty"] = []
    return found, total


def bucketTotal(work):
    bucketDir = os.path.join(work, "buckets")
    if not os.path.isdir(bucketDir):
        return 0
    return sum(os.path.getsize(os.path.join(bucketDir, f)) for f in os.listdir(bucketDir)
               if f.endswith(".bin")) // RECORD_BYTES


def solve(pair, client, curve, work, extra):
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
    proc = subprocess.run(cmd, capture_output=True, text=True)
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
    ap.add_argument("--dp-weight", type=int, default=None,
                    help="the campaign's cutoff; the rewalk must stop at the same points the workers reported")
    ap.add_argument("--client-arg", action="append", default=[],
                    help="extra client argument for the solve step (repeatable)")
    args = ap.parse_args()
    if args.dp_weight is not None:
        args.client_arg = ["--dp-weight", str(args.dp_weight)] + args.client_arg

    os.makedirs(args.work, exist_ok=True)
    statePath = os.path.join(args.work, "state.json")
    state = loadState(statePath)
    if not os.path.exists(statePath):
        if args.buckets & (args.buckets - 1):
            sys.exit("--buckets must be a power of two")
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
    added = ingest(state, root, args.work)
    saveState(statePath, state)
    log("ingested %d new records in %.1f s (%d buckets to re-sort)"
        % (added, time.time() - t0, len(state["dirty"])))
    t1 = time.time()
    found, sorted_ = detect(state, args.work)
    known = {(c["seedA"], c["seedB"]) for c in state["collisions"]}
    new = [c for c in found if (c["seedA"], c["seedB"]) not in known]
    state["collisions"].extend(new)
    saveState(statePath, state)
    # Every recorded pair without a solve result is still pending: a
    # detect-only pass records pairs, a later pass solves them.
    pending = [c for c in state["collisions"] if c.get("result") is None]
    log("sorted %d records in %.1f s; corpus %d records; %d collision(s), %d new, %d unsolved"
        % (sorted_, time.time() - t1, bucketTotal(args.work), len(state["collisions"]), len(new), len(pending)))
    summary = {"added": added, "corpus": bucketTotal(args.work),
               "collisions": len(state["collisions"]), "new": new, "pending": len(pending), "solution": None}
    if args.detect_only or not pending:
        print(json.dumps(summary, indent=1))
        return 0

    for pair in pending:
        res = solve(pair, args.client, args.curve, args.work, args.client_arg)
        pair["result"] = {"k": res["k"], "verified": res["verified"], "returncode": res["returncode"],
                          "tail": res["tail"]}
        saveState(statePath, state)
        log("seeds %s / %s: k = %s, verified %s" % (pair["seedA"], pair["seedB"], res["k"], res["verified"]))
        if res["k"] and res["verified"]:
            res["when"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            state["solved"] = res
            saveState(statePath, state)
            with open(os.path.join(args.work, "solution.json"), "w") as fh:
                json.dump(res, fh, indent=1)
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
