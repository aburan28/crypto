#!/usr/bin/env python3
"""The amendment's admissibility check (RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md §6).

    python3 research/support_local_multipliers_20260924/check_identity.py

Every counter of every rung of every amendment run (postfix/<suite>/<binary>-<arm>/rep*)
must equal the registered run of the same arm (<suite>/<arm>/rep1): the hash fix
changes only how a row is packed, never which rows exist, so any difference
voids the amendment.  Timing fields (*_ns) are excluded; so is nothing else.
Prints one line per suite and arm, then the stage-ladder wall medians of the
four interleaved arms; exits non-zero on any mismatch.
"""
import json, pathlib, statistics, sys

HERE = pathlib.Path(__file__).resolve().parent
SUITES = ["frozen", "chain", "chain-holdout", "chain-holdout-2", "r2-holdout"]
ARMS = ["registered-reference", "registered-candidate", "fixed-reference", "fixed-candidate"]


def counters(row):
    return {k: v for k, v in row.items() if not k.endswith("_ns")}


def main():
    bad = 0
    walls = {}
    for suite in SUITES:
        for arm in ARMS:
            registered = json.loads((HERE / suite / arm.split("-")[1] / "rep1" / "stage.json").read_text())
            reps = sorted((HERE / "postfix" / suite / arm).glob("rep*/stage.json"))
            diffs = 0
            for path in reps:
                doc = json.loads(path.read_text())
                assert len(doc["rows"]) == len(registered["rows"]), path
                for mine, theirs in zip(doc["rows"], registered["rows"]):
                    a, b = counters(mine), counters(theirs)
                    if a != b:
                        diffs += 1
                        keys = sorted(k for k in a.keys() | b.keys() if a.get(k) != b.get(k))
                        print(f"  MISMATCH {path.relative_to(HERE)} {mine['curve']} m={mine['m']}: {keys}")
                walls.setdefault((suite, arm), []).append(sum(r["wall_ns"] for r in doc["rows"]) / 1e9)
            bad += diffs
            print(f"{suite:16s} {arm:22s} reps {len(reps)}  rungs differing from registered: {diffs}")
    print("\nstage-ladder wall, median of three interleaved repetitions (seconds; a practicality note)\n")
    print(f"| suite | {' | '.join(ARMS)} | registered ref/cand | fixed ref/cand | fix, reference arm | fix, candidate arm |")
    print("|:--|" + "--:|" * (len(ARMS) + 4))
    for suite in SUITES:
        m = {arm: statistics.median(walls[(suite, arm)]) for arm in ARMS}
        cells = " | ".join(f"{m[a]:.2f}" for a in ARMS)
        print(
            f"| {suite} | {cells} | {m['registered-reference'] / m['registered-candidate']:.2f}× "
            f"| {m['fixed-reference'] / m['fixed-candidate']:.2f}× "
            f"| {m['registered-reference'] / m['fixed-reference']:.2f}× "
            f"| {m['registered-candidate'] / m['fixed-candidate']:.2f}× |"
        )
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
