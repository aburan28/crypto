#!/usr/bin/env python3
"""Check and tabulate research/env_reads_20260924/run.sh.

    python3 research/env_reads_20260924/check.py > research/env_reads_20260924/check.md

1. Identity: every field but timings (*_ns) of every rung of every stage run,
   both arms, equals the merged default's registered run
   (research/support_local_multipliers_20260924/<suite>/candidate/rep1); every
   whole logarithm verifies [k]G = Q and has the same trials, relations and
   oracle word operations on both arms, seed by seed.
2. Wall, a practicality note: stage-ladder suite medians per arm, and the
   whole-log wall ratio main / change as a geometric mean over seeds (median
   of three repetitions per seed) with a 95% paired bootstrap interval.
Exits non-zero on any identity failure.
"""
import json, math, pathlib, random, statistics, sys

HERE = pathlib.Path(__file__).resolve().parent
REGISTERED = HERE.parent / "support_local_multipliers_20260924"
SUITES = ["frozen", "chain", "chain-holdout", "chain-holdout-2", "r2-holdout"]
ARMS = ["main", "change"]
CELLS = {
    "`K_0/2^13`, seeds 201–210, 301–305": ("K0_2^13", list(range(201, 211)) + list(range(301, 306))),
    "`K_0/2^9`, seeds 201–205": ("K0_2^9", list(range(201, 206))),
}


def counters(row):
    return {k: v for k, v in row.items() if not k.endswith("_ns")}


def main():
    bad = 0
    print("## Identity\n")
    walls = {}
    for suite in SUITES:
        registered = json.loads((REGISTERED / suite / "candidate" / "rep1" / "stage.json").read_text())["rows"]
        for arm in ARMS:
            paths = sorted((HERE / "stage" / suite / arm).glob("rep*/stage.json"))
            diffs = 0
            for p in paths:
                rows = json.loads(p.read_text())["rows"]
                diffs += sum(counters(a) != counters(b) for a, b in zip(rows, registered))
                diffs += abs(len(rows) - len(registered))
                walls.setdefault((suite, arm), []).append(sum(r["wall_ns"] for r in rows) / 1e9)
            diffs += 3 - len(paths)  # a missing repetition is a failure, not a pass
            bad += diffs
            print(f"- {suite}, {arm}: {len(paths)} repetitions, rungs differing from the registered default: {diffs}")
    e2e = {}
    for arm in ARMS:
        for cell, seeds in CELLS.values():
            for s in seeds:
                docs = [json.loads(p.read_text()) for p in sorted((HERE / "e2e" / arm).glob(f"rep*/{cell}_seed{s}.json"))]
                ok = all(d["result"]["verified"] and d["result"]["expected"] == d["result"]["recovered"] for d in docs)
                keys = {(d["counts"]["trials"], d["counts"]["relations"], d["counts"]["f4_word_ops"]) for d in docs}
                if not ok or len(keys) != 1 or len(docs) != 3:
                    bad += 1
                    print(f"- MISMATCH e2e {arm} {cell} seed {s}: verified {ok}, counters {keys}")
                    continue
                e2e[(arm, cell, s)] = (keys.pop(), statistics.median(d["elapsed_seconds"] for d in docs))
    if bad:
        sys.exit(1)
    def registered(cell, s):
        d = json.loads((REGISTERED / "e2e" / "candidate" / f"{cell}_seed{s}.json").read_text())["counts"]
        return (d["trials"], d["relations"], d["f4_word_ops"])

    same = all(
        e2e[("main", c, s)][0] == e2e[("change", c, s)][0] == registered(c, s) for c, ss in CELLS.values() for s in ss
    )
    bad += not same
    print(f"- whole logarithms: 20 seeds x 3 repetitions x 2 arms, all verified; counters equal across arms and to the registered default: {same}")
    print("\n## Wall time (practicality note)\n")
    print("| suite | main (s) | change (s) | main / change |\n|:--|--:|--:|--:|")
    for suite in SUITES:
        a, b = statistics.median(walls[(suite, "main")]), statistics.median(walls[(suite, "change")])
        print(f"| {suite} | {a:.2f} | {b:.2f} | {a / b:.3f}× |")
    rng = random.Random(20260924)
    print("\n| whole logarithms | main → change (sum of per-seed medians) | geometric mean [95% paired bootstrap] |\n|:--|:--|:--|")
    for name, (cell, seeds) in CELLS.items():
        logs = [math.log(e2e[("main", cell, s)][1] / e2e[("change", cell, s)][1]) for s in seeds]
        boots = sorted(math.exp(statistics.fmean(rng.choice(logs) for _ in logs)) for _ in range(10_000))
        wa = sum(e2e[("main", cell, s)][1] for s in seeds)
        wb = sum(e2e[("change", cell, s)][1] for s in seeds)
        print(f"| {name} | {wa:.2f} → {wb:.2f} s | {math.exp(statistics.fmean(logs)):.3f}× [{boots[249]:.3f}, {boots[9749]:.3f}] |")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
