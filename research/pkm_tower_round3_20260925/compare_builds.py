#!/usr/bin/env python3
"""The identity check of note section 12.2: the engine with the compact basis
must reproduce round 2 exactly.

Every row of the replay (`replay/*.jsonl`) is matched, by its system and its
degree bound, with the row round 2 or its cross-check wrote for it
(`research/pkm_tower_round2_20260925/runs/`), and every field but the wall
clock (`ms`, `wall_s`) must be equal. A row without a counterpart, or a
counterpart without a replay, is reported. Every step of every trace that both
builds printed (`tower step k: ...` lines in the `.log` files) must also be
equal on every field round 2 printed, the times aside.

    python3 compare_builds.py            # exits non-zero on any difference
"""

import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND2 = os.path.join(HERE, "..", "pkm_tower_round2_20260925", "runs")
REPLAY = os.path.join(HERE, "replay")

# Wall clock, which the machine decides.
TIMING = {"ms", "wall_s"}
# Round-2 runs of sizes the replay leaves out (note section 12.2).
NOT_REPLAYED = {
    ("M4-kummer-m4-p1", 16),  # ran out of memory in round 2: round 3's cell
    ("M3-kummer-m3-p1", 18),  # stopped for size in round 2
}


def key(r):
    """The system and the degree bound it ran under."""
    return (
        r["p"], r["kind"], r["m"], r["control"], r["t"], r["g"], r["target"],
        r["target_index"], r["x_r"], r["curve"]["a"], r["curve"]["b"],
        json.dumps(r["tower"], sort_keys=True), r.get("max_degree_bound"),
        r.get("staircase_at_stop") is not None,
    )


def rows(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                r = json.loads(line)
                if "summary" not in r and "N" in r:
                    out.append(r)
    return out


def trace(path):
    """`tower step` lines, keyed by their order within each system, without
    the times and without the fields only the new build prints."""
    systems, cur = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("tower step "):
                body = line.split(":", 1)[1]
                # Round 2's fields end at "pairs left N"; the times follow.
                m = re.match(r"(.*?pairs left \d+)", body)
                cur.append(m.group(1).strip() if m else body.strip())
            elif cur:
                systems.append(cur)
                cur = []
    if cur:
        systems.append(cur)
    return systems


def main():
    diffs = 0
    old = {}
    for path in sorted(glob.glob(os.path.join(ROUND2, "*.jsonl"))):
        stem = os.path.basename(path)[: -len(".jsonl")]
        for r in rows(path):
            old.setdefault(key(r), (stem, r))
    compared = identical = 0
    unmatched = []
    replayed = set()
    for path in sorted(glob.glob(os.path.join(REPLAY, "*.jsonl"))):
        stem = os.path.basename(path)[: -len(".jsonl")]
        for r in rows(path):
            k = key(r)
            if k not in old:
                unmatched.append((stem, r["kind"], r["m"], r["control"], r["N"], r["target_index"]))
                continue
            replayed.add(k)
            ostem, o = old[k]
            compared += 1
            fields = (set(r) | set(o)) - TIMING
            bad = sorted(f for f in fields if r.get(f) != o.get(f))
            if bad:
                diffs += 1
                print(f"DIFF {stem} vs {ostem}: {r['kind']} m={r['m']} {r['control']} "
                      f"N={r['N']} target {r['target_index']}: {bad}")
                for f in bad[:6]:
                    print(f"    {f}: replay {r.get(f)!r} / round 2 {o.get(f)!r}")
            else:
                identical += 1
    print(f"rows: {compared} compared, {identical} identical, {compared - identical} different")
    if unmatched:
        diffs += len(unmatched)
        print(f"{len(unmatched)} replay rows have no round-2 counterpart:")
        for u in unmatched[:20]:
            print("   ", u)
    # Round-2 rows the replay should have produced but did not, in the runs
    # that have ended (`progress.txt`).
    ended = set()
    progress = os.path.join(HERE, "progress.txt")
    if os.path.exists(progress):
        with open(progress) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "end":
                    ended.add(parts[2])
    replay_stems = {
        os.path.basename(p)[: -len(".jsonl")] for p in glob.glob(os.path.join(REPLAY, "*.jsonl"))
    }
    running = sorted(replay_stems - ended)
    if running:
        print(f"still running, not checked for missing rows: {running}")
    replay_stems &= ended
    missing = [
        (stem, o["kind"], o["m"], o["control"], o["N"], o["target_index"])
        for k, (stem, o) in old.items()
        if stem in replay_stems and k not in replayed and (stem, o["N"]) not in NOT_REPLAYED
    ]
    if missing:
        diffs += len(missing)
        print(f"{len(missing)} round-2 rows of replayed files were not reproduced:")
        for m in missing[:20]:
            print("   ", m)

    # Traces: every file both builds printed steps into.
    steps = tsame = 0
    for path in sorted(glob.glob(os.path.join(REPLAY, "*.log"))):
        name = os.path.basename(path)
        opath = os.path.join(ROUND2, name)
        if not os.path.exists(opath):
            continue
        new, ref = trace(path), trace(opath)
        if not ref:
            continue
        if len(new) != len(ref):
            diffs += 1
            print(f"TRACE {name}: {len(new)} systems traced, round 2 traced {len(ref)}")
            continue
        for a, b in zip(new, ref):
            if len(a) != len(b):
                diffs += 1
                print(f"TRACE {name}: {len(a)} steps, round 2 had {len(b)}")
                continue
            for i, (x, y) in enumerate(zip(a, b)):
                steps += 1
                if x == y:
                    tsame += 1
                else:
                    diffs += 1
                    print(f"TRACE {name} step {i + 1}:\n    replay  {x}\n    round 2 {y}")
    # Round 2 traced three of these systems under other names: K0's target 0
    # twice while profiling (P0), and M4's N = 12 target 0 (D1). Round 3's
    # cell repeats D2's system, which round 2 stopped after 48 steps: those
    # steps must come out the same.
    others = [
        ("replay", "K0-kummer-m2-p0-N20.log", 0, "P0-profiling-kummer-m2-p0-N20.log", False),
        ("replay", "M4-kummer-m4-p1.log", 2, "D1-kummer-m4-p1-N12-trace.log", False),
        ("runs", "M4b-kummer-m4-p1-N16.log", 0, "D2-kummer-m4-p1-N16-trace.log", True),
    ]
    for where, name, index, ref_name, prefix in others:
        path = os.path.join(HERE, where, name)
        if not os.path.exists(path):
            continue
        new = trace(path)
        if len(new) <= index:
            diffs += 1
            print(f"TRACE {name}: no system {index + 1} traced")
            continue
        a = new[index]
        for k, b in enumerate(trace(os.path.join(ROUND2, ref_name))):
            if len(a) != len(b) and not (prefix and len(a) >= len(b)):
                diffs += 1
                print(f"TRACE {name} against {ref_name} ({k + 1}): {len(a)} steps, round 2 had {len(b)}")
                continue
            for i, (x, y) in enumerate(zip(a, b)):
                steps += 1
                if x == y:
                    tsame += 1
                else:
                    diffs += 1
                    print(f"TRACE {name} against {ref_name} step {i + 1}:\n    new     {x}\n    round 2 {y}")
    print(f"trace steps: {steps} compared, {tsame} identical")
    print("IDENTICAL" if diffs == 0 else f"{diffs} DIFFERENCES")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
