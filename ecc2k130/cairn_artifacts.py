#!/usr/bin/env python3
"""Turn a v2 corpus into cairn claim artifacts.

    python3 cairn_artifacts.py --corpus dp.bin --out artifacts/
    python3 cairn_artifacts.py --corpus dp.bin --limit 64 --stdout

A v2 corpus record is a distinguished point plus its witness: the seed that
started the trail and the eight per-branch step counts (see CAIRN-WITNESS.md).
cairn's orbit-piecework objective pays per orbit for exactly that, as

    {"dps": [{"x": "<orbit>", "seed": "<hex>", "j": [n3, ..., n10]}, ...]}

Two conversions happen here and neither is cosmetic.

The orbit name is **not** this client's `canon`.  Both sides name an orbit by
the least of its members, but in different coordinate orders -- ours is a
permuted type-II ONB, cairn's indexes by gamma^(2^i) -- so the two pick
different members of the same orbit.  `cairn_basis` derives the permutation
between them and this walks the record's canon through it.  Getting that wrong
is silent: the output is a well-formed 131-bit string naming somebody else's
orbit.

The witness is checked before it is shipped.  A record whose counts do not sum
to its `iters` has either overflowed a 32-bit counter or come from a build
with ECC_WITNESS=0, and in both cases the counts are not a witness.  Refusing
here costs nothing; shipping one wastes a claim and, because the unit is the
orbit, takes the payment of whoever reaches it honestly.

What this does NOT do is check that the witness reproduces the endpoint.  That
is one double scalar multiplication over GF(2^131) and it is the checker's
job, on the other side, for its own reasons -- a producer that checked its own
work would be asserting the thing the network exists to verify.

No type hints, camelCase identifiers, no itertools (project convention).
"""

import argparse
import json
import os
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cairn_basis

DP_MAGIC_V2 = b"ECC2KDP2"
DP_HEADER_BYTES = 16
RECORD = struct.Struct("<QQQQQ8I")
JCOUNT = 8


def readCorpus(path):
    """Yield (seed, iters, canon, counts) for every whole v2 record."""
    with open(path, "rb") as fh:
        if fh.read(len(DP_MAGIC_V2)) != DP_MAGIC_V2:
            raise ValueError("%s is not a v2 corpus: no witness in it to submit" % path)
        fh.seek(DP_HEADER_BYTES)
        while True:
            raw = fh.read(RECORD.size)
            if len(raw) < RECORD.size:
                return
            values = RECORD.unpack(raw)
            seed, iters = values[0], values[1]
            canon = values[2] | (values[3] << 64) | (values[4] << 128)
            yield seed, iters, canon, list(values[5:])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--generated", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "generated", "eccF131.h"))
    parser.add_argument("--gamma", default=None,
                        help="the objective's normal-basis generator, polynomial-basis hex")
    parser.add_argument("--batch", type=int, default=64,
                        help="orbits per claim; the objective's max_batch")
    parser.add_argument("--max-weight", dest="maxWeight", type=int, default=34,
                        help="the objective's dp_max_weight; heavier records are not orbits it pays for")
    parser.add_argument("--limit", type=int, default=0, help="stop after this many orbits")
    parser.add_argument("--out", default=None, help="directory to write one JSON file per claim")
    parser.add_argument("--stdout", action="store_true", help="write the claims to stdout")
    args = parser.parse_args()

    zToOnb, sqPerm, _points = cairn_basis.loadGenerated(args.generated)
    gamma = int(args.gamma, 16) if args.gamma else cairn_basis.CAIRN_GAMMA
    try:
        index, perm = cairn_basis.buildPerm(gamma, zToOnb, sqPerm)
    except ValueError as why:
        print("cannot reconcile the bases: %s" % why, file=sys.stderr)
        return 1
    print("basis index %d, %d-bit orbits" % (index, cairn_basis.M), file=sys.stderr)

    if args.out:
        os.makedirs(args.out, exist_ok=True)

    batches = 0
    kept = 0
    witnessless = 0
    heavy = 0
    seen = set()
    duplicates = 0
    pending = []

    def flush():
        nonlocal batches
        if not pending:
            return
        artifact = {"dps": list(pending)}
        text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.out:
            name = os.path.join(args.out, "claim-%05d.json" % batches)
            with open(name, "w") as fh:
                fh.write(text)
        if args.stdout:
            sys.stdout.write(text)
        batches += 1
        del pending[:]

    for seed, iters, canon, counts in readCorpus(args.corpus):
        if sum(counts) != iters:
            witnessless += 1
            continue
        orbit = cairn_basis.leastRotation(cairn_basis.clientToCairn(canon, perm))
        # A corpus also carries the overdue reports the client writes when a
        # lane outruns --max-iters, and those are not distinguished points.
        # The weight is the popcount of the orbit's own coordinates, so this
        # costs nothing and saves shipping a claim the checker will refuse.
        if bin(orbit).count("1") > args.maxWeight:
            heavy += 1
            continue
        name = "%x" % orbit
        # An orbit already claimed is already paid; a second copy verifies and
        # mints nothing, so spending a batch slot on it is pure waste.
        if name in seen:
            duplicates += 1
            continue
        seen.add(name)
        pending.append({"x": name, "seed": "%x" % seed, "j": counts})
        kept += 1
        if len(pending) >= args.batch:
            flush()
        if args.limit and kept >= args.limit:
            break
    flush()

    print("%d orbit(s) in %d claim(s); %d without a usable witness, %d not distinguished, "
          "%d repeated" % (kept, batches, witnessless, heavy, duplicates), file=sys.stderr)
    return 0 if kept else 1


if __name__ == "__main__":
    raise SystemExit(main())
