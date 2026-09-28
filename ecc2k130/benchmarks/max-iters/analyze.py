#!/usr/bin/env python3
"""The analysis PROTOCOL.md fixes: exact step accounting per walk, per-seed
agreement between arms, and the frozen model prediction.

    python3 analyze.py --work DIR [--out DIR]

DIR holds <runId>-<arm>.{bin,ck,log,status} from run.sh.  Writes
results.json and results.txt to --out (default: this directory) and exits
non-zero if any check fails.
"""
import argparse
import hashlib
import json
import math
import os
import re
import struct
import sys

import numpy as np

STEPS = 1024
LAUNCHES = 108
S = STEPS * LAUNCHES
WEIGHT = 44
THETA = 1.447744e-4                 # frozen: even-weight tail at w = 44, m = 131
CAPS = {"ref": 0, "c30": 20800, "c31": 41600, "c32": 83200}
PRODUCTION_CAP = {"ref": None, "c30": 2 ** 30, "c31": 2 ** 31, "c32": 2 ** 32}
RUN_IDS = (30001, 30002, 30003)
HOLDOUTS = (30002, 30003)
THETA_PRODUCTION = 2.0 ** -28.41    # benchmarks/dp-interval
MC_REPLICAS = 400
MC_SEED = 20260925
MC_CHUNK = 50
RARE_TAIL = 0.00135                 # one-sided normal 3-SD tail (amendment 1)
# WALK-CONSTANT.md section 6: iterations per solve on completed trails,
# 2^60.809 x c for the sigma walk's extrapolated constant c in [1.070, 1.082].
RHO_CLASSES_LOG2 = 60.809
SIGMA_C = (1.070, 1.082)

DP_MAGIC = b"ECC2KDP2"
DP_HEADER = 16
DP_RECORD = 72
CK_MAGIC = b"ECC2K130"
CK_HEADER = 40


def idle(t):
    """Steps a walk spends after its trail ends at length t, until the launch
    boundary where the reseed kernel revives it.  Trails start on a boundary."""
    return STEPS * (t // STEPS + 1) - t


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def readCorpus(path):
    blob = open(path, "rb").read()
    if blob[:8] != DP_MAGIC or (len(blob) - DP_HEADER) % DP_RECORD:
        raise ValueError("%s is not a whole v2 corpus" % path)
    n = (len(blob) - DP_HEADER) // DP_RECORD
    rec = np.frombuffer(blob, dtype=np.dtype([
        ("seed", "<u8"), ("iters", "<u8"), ("canon", "<u8", 3), ("counts", "<u4", 8)]),
        count=n, offset=DP_HEADER)
    return rec


def readCheckpoint(path):
    blob = open(path, "rb").read()
    if blob[:8] != CK_MAGIC:
        raise ValueError("%s is not a checkpoint" % path)
    version, m, threads, batch, lanes, runId = struct.unpack_from("<6I", blob, 8)
    iterBase = struct.unpack_from("<Q", blob, 32)[0]
    word = lanes // 8
    words = 2 * threads * batch * m + threads * batch + threads * batch * 8 * 32
    walks = threads * batch * lanes
    if len(blob) != CK_HEADER + words * word + 2 * walks * 8:
        raise ValueError("%s has an unexpected payload size" % path)
    off = CK_HEADER + words * word
    seed = np.frombuffer(blob, dtype="<u8", count=walks, offset=off)
    start = np.frombuffer(blob, dtype="<u8", count=walks, offset=off + walks * 8)
    return dict(version=version, m=m, threads=threads, batch=batch, lanes=lanes,
                runId=runId, iterBase=iterBase, seed=seed, startIter=start)


def parseLog(path):
    text = open(path).read()
    walks = re.search(r"= (\d+) walks, dp weight (\d+), (\d+) steps per launch", text)
    fin = re.search(r"finished: ([\d.]+) M it/s, (\d+) distinguished points "
                    r"\((\d+) verified against the reference, (\d+) dropped\)", text)
    return dict(walks=int(walks.group(1)) if walks else None,
                weight=int(walks.group(2)) if walks else None,
                steps=int(walks.group(3)) if walks else None,
                rate=float(fin.group(1)) if fin else None,
                points=int(fin.group(2)) if fin else None,
                verified=int(fin.group(3)) if fin else None,
                dropped=int(fin.group(4)) if fin else None,
                mismatch="MISMATCH" in text)


def poissonTail(n, lam):
    """P(X >= n) for X ~ Poisson(lam)."""
    term, below = math.exp(-lam), 0.0
    for j in range(n):
        below += term
        term *= lam / (j + 1)
    return max(0.0, 1.0 - below)


def steadyState(cap, theta):
    """Share of trail steps discarded, and share of trails cut, in the long run."""
    log = math.log1p(-theta)
    ac = math.exp(cap * log)
    q = ac * (1.0 - theta)
    kept = (1.0 - theta) / theta * (1.0 - ac * (1.0 + cap * theta))
    return cap * q / (kept + cap * q), q


def simulate(cap, theta, walks, rng):
    """D-hat and q-hat of one run of the exact process, for MC_CHUNK replicas."""
    n = MC_CHUNK * walks
    start = np.zeros(n, dtype=np.int64)
    ended = np.zeros(n, dtype=np.int64)
    cut = np.zeros(n, dtype=np.int64)
    kept = np.zeros(n, dtype=np.int64)
    live = np.arange(n)
    while live.size:
        t = rng.geometric(theta, size=live.size).astype(np.int64) - 1
        isCut = t > cap if cap else np.zeros(live.size, dtype=bool)
        length = np.where(isCut, cap, t)
        event = start[live] + length
        done = event <= S - 1
        idx = live[done]
        ended[idx] += 1
        cut[idx] += isCut[done]
        kept[idx] += np.where(isCut[done], 0, t[done])
        nxt = STEPS * (event[done] // STEPS + 1)
        start[idx] = nxt
        live = idx[nxt < S]
    ended = ended.reshape(MC_CHUNK, walks).sum(axis=1)
    cut = cut.reshape(MC_CHUNK, walks).sum(axis=1)
    kept = kept.reshape(MC_CHUNK, walks).sum(axis=1)
    dropped = cut * cap
    return dropped / (kept + dropped), cut / ended, cut


def analyzeRun(work, runId, arm, failures):
    base = os.path.join(work, "%d-%s" % (runId, arm))
    cap = CAPS[arm]
    status = open(base + ".status").read().split()
    exitCode, seconds = int(status[0]), float(status[1])
    log = parseLog(base + ".log")
    rec = readCorpus(base + ".bin")
    ck = readCheckpoint(base + ".ck")

    def fail(msg):
        failures.append("%d %s: %s" % (runId, arm, msg))

    if exitCode != 0:
        fail("exit %d" % exitCode)
    if log["mismatch"] or log["verified"] != 4 or log["dropped"] != 0:
        fail("verify/dropped: %r" % log)
    if log["points"] != len(rec):
        fail("log reports %s points, corpus holds %d" % (log["points"], len(rec)))
    if (ck["runId"], ck["threads"], ck["iterBase"]) != (runId, 1, S):
        fail("checkpoint header %r" % {k: ck[k] for k in ("runId", "threads", "iterBase")})
    walks = ck["threads"] * ck["batch"] * ck["lanes"]
    if log["walks"] != walks:
        fail("log names %s walks, checkpoint %d" % (log["walks"], walks))

    seeds = rec["seed"].astype(np.uint64)
    iters = rec["iters"].astype(np.int64)
    if np.any((seeds >> np.uint64(48)) != runId):
        fail("a record carries another run id")
    walkIdx = ((seeds >> np.uint64(16)) & np.uint64(0xFFFFFFFF)).astype(np.int64)
    k = (seeds & np.uint64(0xFFFF)).astype(np.int64)
    if len(np.unique(seeds)) != len(seeds):
        fail("a seed is reported twice")
    if np.any(rec["counts"].astype(np.int64).sum(axis=1) != iters):
        fail("witness counts do not sum to the trail length")
    if cap and np.any(iters > cap):
        fail("a reported trail is longer than the cap")

    ckSeed = ck["seed"].astype(np.uint64)
    lane = np.arange(walks, dtype=np.uint64)
    if np.any(ckSeed >> np.uint64(16) != (np.uint64(runId) << np.uint64(32)) | lane):
        fail("checkpoint seeds are not this run's walk indices")
    endedPerWalk = (ckSeed & np.uint64(0xFFFF)).astype(np.int64)
    if np.any(k >= endedPerWalk[walkIdx]):
        fail("a record's restart index is past its walk's counter")
    reportedPerWalk = np.bincount(walkIdx, minlength=walks)
    cutPerWalk = endedPerWalk - reportedPerWalk
    if np.any(cutPerWalk < 0) or (cap == 0 and np.any(cutPerWalk)):
        fail("cut count negative, or cuts without a cap")
    trailSteps = np.bincount(walkIdx, weights=iters + idle(iters),
                             minlength=walks).astype(np.int64)
    start = ck["startIter"].astype(np.int64)
    expect = trailSteps + cutPerWalk * ((cap + idle(cap)) if cap else 0)
    bad = int(np.count_nonzero(expect != start))
    if bad:
        fail("accounting identity fails on %d walks" % bad)

    kept = int(iters.sum())
    ncut = int(cutPerWalk.sum())
    dropped = ncut * cap
    idleSteps = int(idle(iters).sum()) + ncut * (idle(cap) if cap else 0)
    inflight = int((S - start).sum())
    if kept + dropped + idleSteps + inflight != walks * S:
        fail("steps do not add up to walks x S")
    ended = int(endedPerWalk.sum())
    return dict(
        runId=runId, arm=arm, cap=cap, exit=exitCode, seconds=seconds, rate=log["rate"],
        walks=walks, lanes=ck["lanes"], ended=ended, reported=len(rec), cut=ncut,
        keptSteps=kept, droppedSteps=dropped, idleSteps=idleSteps, inflightSteps=inflight,
        D=dropped / (kept + dropped), q=ncut / ended, loss=(kept + dropped) / kept,
        censoredChecks=int((iters + 1).sum() + inflight), verified=log["verified"],
        corpus={"bytes": os.path.getsize(base + ".bin"), "sha256": sha256(base + ".bin")},
        checkpoint={"bytes": os.path.getsize(base + ".ck"), "sha256": sha256(base + ".ck")},
        _seeds=seeds, _iters=iters, _canon=rec["canon"], _ended=endedPerWalk)


def crossCheck(runs, runId, failures):
    """Per seed, the arms of one run id agree (PROTOCOL.md check 3)."""
    ref = runs[(runId, "ref")]
    order = np.argsort(ref["_seeds"])
    refSeeds, refIters = ref["_seeds"][order], ref["_iters"][order]
    refCanon = ref["_canon"][order]
    out = {}
    for arm in ("c30", "c31", "c32"):
        r = runs[(runId, arm)]
        cap = CAPS[arm]
        pos = np.searchsorted(refSeeds, r["_seeds"])
        pos = np.minimum(pos, len(refSeeds) - 1)
        both = refSeeds[pos] == r["_seeds"]
        sameLen = np.all(refIters[pos[both]] == r["_iters"][both])
        sameOrbit = np.all(refCanon[pos[both]] == r["_canon"][both])
        if not (sameLen and sameOrbit):
            failures.append("%d %s: a seed reported by both arms differs" % (runId, arm))
        # A reference trail longer than the cap is never reported under it; one
        # within the cap is reported whenever the capped arm finished its seed.
        refWalk = ((refSeeds >> np.uint64(16)) & np.uint64(0xFFFFFFFF)).astype(np.int64)
        refK = (refSeeds & np.uint64(0xFFFF)).astype(np.int64)
        finished = refK < r["_ended"][refWalk]
        inArm = np.isin(refSeeds, r["_seeds"])
        long = refIters > cap
        wrongReport = int(np.count_nonzero(long & inArm))
        wrongCut = int(np.count_nonzero(~long & finished & ~inArm))
        if wrongReport or wrongCut:
            failures.append("%d %s: %d long trails reported, %d short trails cut"
                            % (runId, arm, wrongReport, wrongCut))
        out[arm] = dict(sharedSeeds=int(both.sum()), identical=bool(sameLen and sameOrbit),
                        refLongFinished=int(np.count_nonzero(long & finished)),
                        refShortFinished=int(np.count_nonzero(~long & finished)),
                        longReported=wrongReport, shortCut=wrongCut)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    args = ap.parse_args()

    failures = []
    runs = {}
    for runId in RUN_IDS:
        for arm in CAPS:
            runs[(runId, arm)] = analyzeRun(args.work, runId, arm, failures)
    cross = {runId: crossCheck(runs, runId, failures) for runId in RUN_IDS}

    refs = [runs[(r, "ref")] for r in RUN_IDS]
    thetaHat = sum(r["reported"] for r in refs) / sum(r["censoredChecks"] for r in refs)
    thetaSe = math.sqrt(THETA * (1 - THETA) / sum(r["censoredChecks"] for r in refs))

    walks = runs[(RUN_IDS[0], "ref")]["walks"]
    rng = np.random.default_rng(MC_SEED)
    model = {}
    for arm in ("c30", "c31", "c32"):
        cap = CAPS[arm]
        ds, qs, cs = [], [], []
        for _ in range(MC_REPLICAS // MC_CHUNK):
            d, q, c = simulate(cap, THETA, walks, rng)
            ds.extend(d.tolist())
            qs.extend(q.tolist())
            cs.extend(c.tolist())
        dss, qss = steadyState(cap, THETA)
        model[arm] = dict(Dmean=float(np.mean(ds)), Dsd=float(np.std(ds, ddof=1)),
                          qmean=float(np.mean(qs)), qsd=float(np.std(qs, ddof=1)),
                          cutMean=float(np.mean(cs)), Dsteady=dss, qsteady=qss)

    verdict = {}
    for arm in ("c30", "c31", "c32"):
        m = model[arm]
        for runId in RUN_IDS:
            r = runs[(runId, arm)]
            zD = (r["D"] - m["Dmean"]) / m["Dsd"] if m["Dsd"] > 0 else 0.0
            zq = (r["q"] - m["qmean"]) / m["qsd"] if m["qsd"] > 0 else 0.0
            r["zD"], r["zq"] = zD, zq
            if arm == "c32":
                # PROTOCOL.md amendment 1: under one cut per run is a Poisson
                # count, judged by its exact tail rather than a z-score.
                r["tail"] = poissonTail(r["cut"], m["cutMean"])
                ok = r["tail"] >= RARE_TAIL and r["D"] < 1e-3
            else:
                ok = abs(zD) <= 3 and abs(zq) <= 3
            verdict["%d %s" % (runId, arm)] = ok
            if not ok:
                failures.append("%d %s: outside the model (zD %.2f, zq %.2f, tail %s)"
                                % (runId, arm, zD, zq, r.get("tail")))

    production = {}
    for arm, cap in PRODUCTION_CAP.items():
        if cap is None:
            continue
        d, q = steadyState(cap, THETA_PRODUCTION)
        loss = 1.0 / (1.0 - d)
        production[arm] = dict(cap=cap, trailLengths=cap * THETA_PRODUCTION, q=q, D=d, loss=loss,
                               workLog2=[RHO_CLASSES_LOG2 + math.log2(c * loss) for c in SIGMA_C])

    public = {
        "protocol": "PROTOCOL.md", "theta": THETA, "thetaHat": thetaHat, "thetaSe": thetaSe,
        "S": S, "steps": STEPS, "weight": WEIGHT, "caps": CAPS, "runIds": RUN_IDS,
        "holdouts": HOLDOUTS, "mc": {"replicas": MC_REPLICAS, "seed": MC_SEED},
        "runs": [{k: v for k, v in r.items() if not k.startswith("_")} for r in runs.values()],
        "cross": {str(k): v for k, v in cross.items()},
        "model": model, "verdict": verdict, "production": production,
        "failures": failures, "passed": not failures,
    }
    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump(public, fh, indent=1, sort_keys=True)
        fh.write("\n")
    text = report(public, runs)
    with open(os.path.join(args.out, "results.txt"), "w") as fh:
        fh.write(text)
    sys.stdout.write(text)
    return 0 if not failures else 1


def report(p, runs):
    out = []
    w = out.append
    w("scale model: curve 131, weight %d, %d walks, S = %d steps per walk, %d-step launches"
      % (p["weight"], runs[(RUN_IDS[0], "ref")]["walks"], S, STEPS))
    w("theta frozen %.6e, measured on the reference arms %.6e (ratio %.4f, %.1f standard errors)"
      % (THETA, p["thetaHat"], p["thetaHat"] / THETA, (p["thetaHat"] - THETA) / p["thetaSe"]))
    w("")
    w("%-4s %-6s %6s %6s %8s %7s %7s %9s %9s %9s %6s %6s %8s %s"
      % ("arm", "run", "cap", "c/mu", "ended", "cut", "q%", "D%", "model D%", "sd", "zD", "zq",
         "loss", "checks"))
    for arm in CAPS:
        for runId in RUN_IDS:
            r = runs[(runId, arm)]
            m = p["model"].get(arm)
            tag = "dev" if runId == RUN_IDS[0] else "hold"
            capText = "inf" if not r["cap"] else str(r["cap"])
            ratio = "inf" if not r["cap"] else "%.3f" % (r["cap"] * THETA)
            if m:
                w("%-4s %-6s %6s %6s %8d %7d %7.3f %9.4f %9.4f %9.4f %6.2f %6.2f %8.5f %s"
                  % (arm, tag, capText, ratio, r["ended"], r["cut"], 100 * r["q"], 100 * r["D"],
                     100 * m["Dmean"], 100 * m["Dsd"], r["zD"], r["zq"], r["loss"],
                     "ok" if p["verdict"]["%d %s" % (runId, arm)] else "FAIL"))
            else:
                w("%-4s %-6s %6s %6s %8d %7d %7.3f %9.4f %9s %9s %6s %6s %8.5f %s"
                  % (arm, tag, capText, ratio, r["ended"], r["cut"], 100 * r["q"], 100 * r["D"],
                     "-", "-", "-", "-", r["loss"], "reference"))
    w("")
    w("steady state at the frozen theta: " + ", ".join(
        "%s D = %.4f%% q = %.4f%%" % (a, 100 * m["Dsteady"], 100 * m["qsteady"])
        for a, m in p["model"].items()))
    w("")
    w("production projection (model row: theta = 2^-28.41, guard overshoot and launch idle omitted)")
    for arm, pr in p["production"].items():
        w("  maxIters 2^%d = %.2f trail lengths: cuts %.4f%% of trails, discards %.5f%% of steps, "
          "loss x%.5f, sigma work 2^%.2f-%.2f"
          % (round(math.log2(pr["cap"])), pr["trailLengths"], 100 * pr["q"], 100 * pr["D"],
             pr["loss"], pr["workLog2"][0], pr["workLog2"][1]))
    w("")
    for runId, arms in p["cross"].items():
        w("per-seed agreement, run %s: " % runId + "; ".join(
            "%s %d shared seeds identical=%s, %d long / %d short reference trails finished, "
            "%d long reported, %d short cut" % (a, v["sharedSeeds"], v["identical"],
                                                 v["refLongFinished"], v["refShortFinished"],
                                                 v["longReported"], v["shortCut"])
            for a, v in arms.items()))
    w("")
    w("PASSED" if p["passed"] else "FAILED:\n  " + "\n  ".join(p["failures"]))
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.exit(main())
