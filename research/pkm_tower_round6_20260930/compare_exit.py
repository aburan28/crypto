#!/usr/bin/env python3
"""The identity check of note section 15.2: the engine with the full-rank exit
against the rows and traces rounds 2 and 3 committed for the same systems.

    python3 research/pkm_tower_round6_20260930/compare_exit.py

- `off/` (the exit off) must reproduce every field of every row but the wall
  clock (`ms`, `wall_s`): the multiply-adds too.
- `replay/` (the exit on) must reproduce every field but the wall clock, the
  multiply-adds and `max_residual_rows`, and those two may only fall.
- Every step of every trace both builds printed must be equal on every field
  the old build printed, the times aside, but for the residue count, which
  may only fall.

Fields only the new build writes (`full_rank_exit`, `rows_skipped_full_rank`)
are ignored. A replay row without a counterpart, or a counterpart of a
replayed file without a replay, is reported. The multiply-adds saved are
listed per system. Exits non-zero on any difference.
"""

import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = [
    os.path.join(HERE, "..", "pkm_tower_round2_20260925", "runs"),
    os.path.join(HERE, "..", "pkm_tower_round3_20260925", "runs"),
]
TIMING = {"ms", "wall_s"}
NEW = {"full_rank_exit", "rows_skipped_full_rank"}
MAY_FALL = {"muladds", "max_residual_rows"}
# Round-2 sizes the replay leaves out, as round 3's did (note section 12.2).
NOT_REPLAYED = {("M4-kummer-m4-p1", 16), ("M3-kummer-m3-p1", 18)}


def key(r):
    """The system and the degree bound it ran under (round 3's rule)."""
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


RESIDUE = re.compile(r"residue (\d+) x ")


def trace(path):
    """`tower step` lines per system: the fields up to `pairs left N` with the
    residue count taken out, and the residue count."""
    systems, cur = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("tower step "):
                body = line.split(":", 1)[1]
                m = re.match(r"(.*?pairs left \d+)", body)
                fields = (m.group(1) if m else body).strip()
                r = RESIDUE.search(fields)
                residue = int(r.group(1)) if r else None
                cur.append((RESIDUE.sub("residue * x ", fields), residue, line))
            elif cur:
                systems.append(cur)
                cur = []
    if cur:
        systems.append(cur)
    return systems


def ended():
    done = set()
    path = os.path.join(HERE, "progress.txt")
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "end":
                    done.add(parts[2])
    return done


def main():
    diffs = 0
    ref = {}
    for d in REFS:
        for path in sorted(glob.glob(os.path.join(d, "*.jsonl"))):
            stem = os.path.basename(path)[: -len(".jsonl")]
            for r in rows(path):
                ref.setdefault(key(r), (stem, r))
    done = ended()
    savings = []
    total = 0
    for where, may_fall in (("off", set()), ("replay", MAY_FALL)):
        compared = identical = 0
        seen, stems = set(), set()
        for path in sorted(glob.glob(os.path.join(HERE, where, "*.jsonl"))):
            stem = os.path.basename(path)[: -len(".jsonl")]
            if f"{where}/{stem}" not in done:
                print(f"{where}/{stem}: still running, not compared")
                continue
            stems.add(stem)
            for r in rows(path):
                k = key(r)
                if k not in ref:
                    diffs += 1
                    print(f"{where}/{stem}: no reference row for {r['kind']} m={r['m']} "
                          f"{r['control']} N={r['N']} target {r['target_index']}")
                    continue
                seen.add(k)
                ostem, o = ref[k]
                compared += 1
                fields = (set(r) | set(o)) - TIMING - NEW
                bad = [f for f in sorted(fields) if f not in may_fall and r.get(f) != o.get(f)]
                bad += [f for f in sorted(may_fall) if (r.get(f) or 0) > (o.get(f) or 0)]
                if bad:
                    diffs += 1
                    print(f"DIFF {where}/{stem} vs {ostem}: {r['kind']} m={r['m']} {r['control']} "
                          f"N={r['N']} target {r['target_index']}: {bad}")
                    for f in bad[:6]:
                        print(f"    {f}: new {r.get(f)!r} / committed {o.get(f)!r}")
                    continue
                identical += 1
                if where == "replay" and r["muladds"] < o["muladds"]:
                    savings.append((stem, r["kind"], r["m"], r["control"], r["N"], r["target_index"],
                                    o["muladds"], r["muladds"], r.get("rows_skipped_full_rank")))
        total += compared
        print(f"{where}: {compared} rows compared, {identical} as required")
        missing = [
            (stem, o["kind"], o["m"], o["control"], o["N"], o["target_index"])
            for k, (stem, o) in ref.items()
            if stem in stems and k not in seen and (stem, o["N"]) not in NOT_REPLAYED
        ]
        if missing:
            diffs += len(missing)
            print(f"{len(missing)} committed rows of {where} files were not reproduced:")
            for m in missing[:20]:
                print("   ", m)

        # Traces.
        steps = same = exits = 0
        for path in sorted(glob.glob(os.path.join(HERE, where, "*.log"))):
            name = os.path.basename(path)
            if f"{where}/{name[: -len('.log')]}" not in done:
                continue
            opath = next((os.path.join(d, name) for d in REFS if os.path.exists(os.path.join(d, name))), None)
            if opath is None:
                continue
            new, old = trace(path), trace(opath)
            if not old:
                continue
            if len(new) != len(old):
                diffs += 1
                print(f"TRACE {where}/{name}: {len(new)} systems traced, committed {len(old)}")
                continue
            for a, b in zip(new, old):
                if len(a) != len(b):
                    diffs += 1
                    print(f"TRACE {where}/{name}: {len(a)} steps, committed {len(b)}")
                    continue
                for i, ((x, rx, line), (y, ry, _)) in enumerate(zip(a, b)):
                    steps += 1
                    ok = x == y and (rx == ry if where == "off" else (rx or 0) <= (ry or 0))
                    if ok:
                        same += 1
                    else:
                        diffs += 1
                        print(f"TRACE {where}/{name} step {i + 1}:\n    new       {x} (residue {rx})"
                              f"\n    committed {y} (residue {ry})")
                    exits += int("skipped)" in line and "(0 skipped)" not in line)
        print(f"{where} trace steps: {steps} compared, {same} as required, {exits} with rows skipped")

    print("\nMultiply-adds saved by the exit (replay against the committed rows):\n")
    print("| file | kind | m | control | N | target | committed | replay | saved | S-rows skipped |")
    print("|:--|:--|--:|:--|--:|--:|--:|--:|--:|--:|")
    for stem, kind, m, control, n, t, old, new, skipped in savings:
        print(f"| {stem} | {kind} | {m} | {control} | {n} | {t} | {old:.4g} | {new:.4g} | "
              f"{100 * (old - new) / old:.1f}% | {skipped} |")
    if total == 0 and diffs == 0:
        print("\nNo round-6 rows yet: nothing compared.")
    else:
        print("\nIDENTICAL, as required" if diffs == 0 else f"\n{diffs} DIFFERENCES")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
