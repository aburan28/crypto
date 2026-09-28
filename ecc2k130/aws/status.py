#!/usr/bin/env python3
"""Campaign dashboard from the slot table: who is walking, how fast, how far.

    status.py                 one snapshot
    status.py --watch 60      refresh every 60 s
    status.py --json          machine-readable

Progress is measured against 2^60.9 expected iterations (Bailey et al.),
which is an expectation, not a deadline: the probability that a collision has
not happened after c times the expected work is about exp(-pi c^2 / 4)
(17% at c = 1.5, 4% at c = 2).  Iterations counted here are the sum of every
slot's checkpointed progress times that slot's own walk count, so they
survive worker restarts; a slot's live-but-unuploaded stretch is not counted.

A checkpoint stores the per-walk iteration base, so turning it into group
operations needs the number of walks the slot is actually running.  That is
not a campaign constant: Ada slots omit --threads and let the client size the
grid (worker.usesCampaignWorkers), so they carry far fewer walks than the
385,024 x 16 Blackwell preset and their iteration base climbs proportionally
faster.  Multiplying every slot by the preset inflated such a slot's
contribution by the ratio of the two grids -- tens of times over, with no
change in the points it produced.  Workers now report `walks`; slots that do
not (an old worker) are counted at --walks and listed as assumed, so the
headline says how much of itself it is guessing.

No type hints, camelCase identifiers (project convention).
"""

import argparse
import json
import os
import subprocess
import sys
import time

EXPECTED_ITERS = 2.0 ** 60.9
LEASE_GRACE = 0


def fromDv(v):
    if "N" in v:
        n = v["N"]
        return float(n) if "." in n or "e" in n else int(n)
    if "S" in v:
        return v["S"]
    if "BOOL" in v:
        return v["BOOL"]
    return None


def scan(table, local=None, bucket=None):
    if local:
        with open(os.path.join(local, "slots.json")) as fh:
            return [dict(v, slot=int(k)) for k, v in json.load(fh).items()]
    if table:
        r = subprocess.run(["aws", "dynamodb", "scan", "--table-name", table, "--output", "json"],
                           capture_output=True, text=True, check=True)
        return [{k: fromDv(v) for k, v in it.items()} for it in json.loads(r.stdout).get("Items", [])]
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from worker import S3Slots
    return S3Slots(bucket).scan()


def humanTime(sec):
    if sec < 3600:
        return "%dm" % (sec // 60)
    if sec < 86400:
        return "%.1fh" % (sec / 3600)
    return "%.1fd" % (sec / 86400)


def slotWalks(it, fallback):
    """The walks this slot is running: what it reported, else the preset."""
    walks = int(it.get("walks") or 0)
    return walks if walks > 0 else fallback


def slotIters(it, fallback):
    return max(0, int(it.get("ckptIter") or 0)) * slotWalks(it, fallback)


def report(items, walksPerSlot, asJson):
    now = time.time()
    alive = [it for it in items if int(it.get("leaseUntil") or 0) >= now]
    rate = sum(float(it.get("rate") or 0) for it in alive)
    ckptIters = sum(slotIters(it, walksPerSlot) for it in items)
    assumed = [it for it in items
               if int(it.get("walks") or 0) <= 0 and int(it.get("ckptIter") or 0) > 0]
    dps = sum(int(it.get("dpUploaded") or 0) for it in items)
    solved = [it for it in items if it.get("state") == "solved"]
    frac = ckptIters / EXPECTED_ITERS
    remaining = (EXPECTED_ITERS - ckptIters) / rate if rate > 0 else float("inf")
    summary = {"slots": len(items), "alive": len(alive),
               "retired": sum(1 for it in items if it.get("state") == "retired"),
               "errors": sum(1 for it in items if it.get("state") == "error"),
               "rateBps": rate, "checkpointedIterations": ckptIters, "fractionOfExpected": frac,
               "dpUploaded": dps, "etaSecondsAtCurrentRate": remaining, "solved": solved,
               "walksAssumedSlots": len(assumed), "assumedWalksPerSlot": walksPerSlot,
               "checkpointedIterationsAssumed": sum(slotIters(it, walksPerSlot) for it in assumed)}
    if asJson:
        print(json.dumps(summary, indent=1, default=str))
        return
    print("slots %d, alive %d, retired %d, errors %d" % (len(items), len(alive), summary["retired"], summary["errors"]))
    print("aggregate %.3f B it/s   checkpointed %.4g iterations = %.3f%% of 2^60.9   %d points uploaded"
          % (rate / 1e9, ckptIters, 100 * frac, dps))
    if assumed:
        print("  %d slot(s) do not report their walk count; counted at the %d-walk preset "
              "(%.4g of the iterations above are that guess)"
              % (len(assumed), walksPerSlot, summary["checkpointedIterationsAssumed"]))
    if rate > 0:
        print("at this rate the expected remaining work takes %s (17%% chance it needs 1.5x, 4%% chance 2x)"
              % humanTime(remaining))
    if solved:
        print("SOLVED: %s" % solved[0].get("solution"))
    print("%5s %-28s %-24s %7s %10s %12s %10s %8s" % (
        "slot", "owner", "gpu", "B it/s", "dp", "ckpt iter", "walks", "lease"))
    for it in sorted(items, key=lambda x: x["slot"]):
        lease = int(it.get("leaseUntil") or 0) - now
        state = it.get("state", "?")
        status = "%ds" % lease if lease > 0 else ("expired" if state == "active" else state)
        walks = int(it.get("walks") or 0)
        print("%5d %-28s %-24s %7.3f %10d %12d %10s %8s" % (
            it["slot"], str(it.get("owner", "-"))[:28], str(it.get("gpuName", "-"))[:24],
            float(it.get("rate") or 0) / 1e9, int(it.get("dpUploaded") or 0),
            int(it.get("ckptIter") or 0), walks if walks > 0 else "%d?" % walksPerSlot, status))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bucket", default=os.environ.get("ECC_BUCKET"), help="campaign bucket (slots/ prefix)")
    ap.add_argument("--table", default=os.environ.get("ECC_TABLE", ""), help="DynamoDB table, if slots live there")
    ap.add_argument("--local", default=os.environ.get("ECC_LOCAL_STORE"), help="rehearsal store directory")
    ap.add_argument("--walks", type=int, default=385024 * 16,
                    help="walks assumed for a slot that does not report its own (workers x batch)")
    ap.add_argument("--watch", type=float, default=0)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    if not (args.local or args.table or args.bucket):
        sys.exit("give --bucket (or ECC_BUCKET), --table or --local")
    while True:
        try:
            report(scan(args.table, args.local, args.bucket), args.walks, args.json)
        except Exception as e:
            print("status failed: %s" % e, file=sys.stderr)
        if not args.watch:
            return 0
        time.sleep(args.watch)
        print()


if __name__ == "__main__":
    sys.exit(main())
