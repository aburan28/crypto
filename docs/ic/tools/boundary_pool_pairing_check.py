#!/usr/bin/env python3
"""Check that drawing the restart pool lazily moved the cost and nothing else.

usage: boundary_pool_pairing_check.py <eager.json> <lazy.json> [--out <json>]

The note's §13 claims something stronger than "the candidate is faster":
it claims the two arms walk *bit-identical trajectories* and differ in
exactly one quantity, the offsets each row paid for.  That is checkable
rather than arguable, and this checks it.

For each row paired by `(regime, instance, variant, repeat)`:

  * every native counter and group-operation count must be identical,
    excepting only the three that report the pool itself
    (`walk_jumps`, `walk_pool_offsets`, `walk_pool_ops`) and the
    `relations.group_ops` fields the offsets are charged to;
  * the saving must account for itself exactly — the drop in the
    relation phase's group-addition equivalents must *equal* the drop in
    `walk_pool_ops`, and the drop in scalar multiplications must be
    twice the drop in offsets;
  * no row's `S` may rise.

It then reports where the saving actually lands, in `S` and as a ratio,
because §13.1 predicts it is confined to the rungs below `2^20`.

Wall-clock fields are excluded: they never priced anything.
"""
import json, math, sys
from collections import OrderedDict

PHASES = ("factor_base", "relations", "linear_algebra", "verify")
# The counters that report the pool, and the ledger fields it is charged
# to.  Everything else must match to the operation.
POOL_KEYS = {
    "relations.native.walk_jumps",
    "relations.native.walk_pool_offsets",
    "relations.native.walk_pool_ops",
    "relations.group_ops.adds",
    "relations.group_ops.doubles",
    "relations.group_ops.scalar_mults",
}


def rows_of(run):
    out = {}
    for inst in run["ledger"]["instances"]:
        seen = OrderedDict()
        for v in inst["variants"]:
            seen.setdefault(v["name"], []).append(v)
        out[(inst["regime"], inst["curve"]["name"])] = (inst, seen)
    return out


def counters(v):
    flat = {"trials": v["trials"], "relations_found": v["relations_found"],
            "independent": v["independent"], "dependent": v["dependent"],
            "rank": v["rank"], "recovered": v.get("recovered"),
            "verified": v["verified"], "signed_points": v["signed_points"],
            "columns": v["columns"], "summands": v["summands"]}
    for p in PHASES:
        if p not in v:
            continue
        for k, x in v[p]["group_ops"].items():
            flat[f"{p}.group_ops.{k}"] = x
        for k, x in v[p]["native"].items():
            flat[f"{p}.native.{k}"] = x
    return flat


def native(v, key, default=0):
    return v["relations"]["native"].get(key, default)


def gae_ops(v):
    """Group-addition equivalents of a phase's own counts: adds + doubles."""
    g = v["relations"]["group_ops"]
    return g["adds"] + g["doubles"]


eager = json.load(open(sys.argv[1]))
lazy = json.load(open(sys.argv[2]))
out_path = next((sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--out"), None)

A, B = rows_of(eager), rows_of(lazy)
identical = 0
trajectory_moved, accounting_broken, risers, residue = [], [], [], []
walk_rows, savers, savings = 0, 0, []
per_row = []
# Instances whose conversion ratios are not in the pinned table price their
# non-addition units at factors measured on the host each run, so two runs
# of such an instance differ by the drift §12 measured even when every
# counter is identical.  That is not this round's doing and is reported
# apart from it.
unpinned = {inst["curve"]["name"]: inst.get("calibration_pinned", {}).get("measured", [])
            for inst in lazy["ledger"]["instances"]
            if inst.get("calibration_pinned", {}).get("measured")}

for key in B:
    if key not in A:
        continue
    (inst_a, rows_a), (inst_b, rows_b) = A[key], B[key]
    log2r = inst_b["log2_r"]
    for name, runs_b in rows_b.items():
        runs_a = rows_a.get(name)
        if not runs_a:
            continue
        for i in range(min(len(runs_a), len(runs_b))):
            va, vb = runs_a[i], runs_b[i]
            ca, cb = counters(va), counters(vb)
            bad = [k for k in set(ca) | set(cb)
                   if k not in POOL_KEYS and ca.get(k) != cb.get(k)]
            if bad:
                trajectory_moved.append((key, name, i, sorted(bad)[:4]))
            else:
                identical += 1
            is_walk = "walk_pool_offsets" in vb["relations"]["native"]
            if is_walk:
                walk_rows += 1
                pool_saving = native(va, "walk_pool_ops") - native(vb, "walk_pool_ops")
                gae_saving = gae_ops(va) - gae_ops(vb)
                mult_saving = (va["relations"]["group_ops"]["scalar_mults"]
                               - vb["relations"]["group_ops"]["scalar_mults"])
                offset_saving = native(va, "walk_pool_offsets") - native(vb, "walk_pool_offsets")
                if gae_saving != pool_saving or mult_saving != 2 * offset_saving:
                    accounting_broken.append(
                        (key, name, i, gae_saving, pool_saving, mult_saving, offset_saving))
                if pool_saving > 0:
                    savers += 1
                    savings.append(pool_saving)
            if va["s"] and vb["s"]:
                ratio = vb["s"] / va["s"]
                delta = va["s"] - vb["s"]
                per_row.append((ratio, delta, log2r, key[0], key[1], name, is_walk))
                if ratio > 1 + 1e-12:
                    (residue if (not bad and key[1] in unpinned) else risers).append(
                        (ratio, key, name, i))

print(f"rows compared: {identical + len(trajectory_moved)}")
print(f"identical in every counter but the pool's: {identical}")
print(f"rows whose trajectory moved: {len(trajectory_moved)}")
for key, name, i, bad in trajectory_moved[:8]:
    print(f"   {key[0]} {key[1]} {name} repeat {i}: {bad}")
print(f"walk rows: {walk_rows}, of which drew fewer offsets: {savers}")
print(f"rows whose saving does not equal the offsets they stopped drawing: {len(accounting_broken)}")
for key, name, i, g, p, m, o in accounting_broken[:8]:
    print(f"   {key[0]} {key[1]} {name} repeat {i}: gae {g} vs pool {p}, mults {m} vs 2x{o}")
print(f"rows whose S rose: {len(risers)}")
for ratio, key, name, i in sorted(risers, key=lambda x: -x[0])[:8]:
    print(f"   {ratio:.4f}  {key[0]} {key[1]} {name} repeat {i}")
if residue:
    print(f"rows that rose on identical counters, all on instances the pinned table "
          f"does not carry: {len(residue)}")
    for inst, units in unpinned.items():
        print(f"   {inst}: measured, not pinned: {', '.join(units)}")

if per_row:
    walk_only = [p for p in per_row if p[6]]
    print()
    print("where the saving lands (walk rows, S before minus S after):")
    for ratio, delta, log2r, regime, inst, name, _ in sorted(walk_only, key=lambda x: -x[1])[:12]:
        print(f"   log2 r {log2r:5.1f}  dS {delta:9.3f}  ratio {ratio:6.4f}  "
              f"{regime:8s} {inst[:24]:24s} {name}")
    inside = [p for p in walk_only if p[2] >= 20]
    outside = [p for p in walk_only if p[2] < 20]
    for label, rows in (("below 2^20", outside), ("at or above 2^20", inside)):
        if rows:
            worst = min(r[0] for r in rows)
            print(f"   {label}: {len(rows)} walk rows, best ratio {worst:.4f}, "
                  f"largest dS {max(r[1] for r in rows):.3f}")

if out_path:
    walk_only = [p for p in per_row if p[6]]
    doc = OrderedDict([
        ("what", "Check that drawing the walk's restart pool on first use moved cost and "
                 "nothing else: every counter of every row compared between the eager and "
                 "the lazy arm, the saving checked against the offsets each row stopped "
                 "drawing, and the size of the saving by instance."),
        ("eager", sys.argv[1]),
        ("lazy", sys.argv[2]),
        ("rows_compared", identical + len(trajectory_moved)),
        ("rows_identical_but_for_the_pool", identical),
        ("rows_whose_trajectory_moved", len(trajectory_moved)),
        ("trajectory_mismatches", [{"regime": k[0], "instance": k[1], "variant": n,
                                    "repeat": i, "fields": b}
                                   for k, n, i, b in trajectory_moved]),
        ("walk_rows", walk_rows),
        ("walk_rows_that_drew_fewer_offsets", savers),
        ("rows_whose_saving_is_not_the_offsets_they_dropped", len(accounting_broken)),
        ("accounting_mismatches", [{"regime": k[0], "instance": k[1], "variant": n, "repeat": i,
                                    "gae_saving": g, "pool_ops_saving": p,
                                    "scalar_mult_saving": m, "offset_saving": o}
                                   for k, n, i, g, p, m, o in accounting_broken]),
        ("rows_whose_s_rose", len(risers)),
        ("rows_that_rose_on_identical_counters_at_unpinned_instances", len(residue)),
        ("instances_the_pinned_table_does_not_carry", unpinned),
        ("largest_savings", [{"log2_r": round(l, 3), "delta_s": round(d, 6),
                              "ratio": round(r, 6), "regime": g, "instance": i, "variant": n}
                             for r, d, l, g, i, n, _ in sorted(walk_only, key=lambda x: -x[1])[:40]]),
        ("best_ratio_below_2_20", round(min([p[0] for p in walk_only if p[2] < 20], default=float("nan")), 6)),
        ("best_ratio_at_or_above_2_20", round(min([p[0] for p in walk_only if p[2] >= 20], default=float("nan")), 6)),
        ("largest_delta_s_at_or_above_2_20", round(max([p[1] for p in walk_only if p[2] >= 20], default=float("nan")), 6)),
    ])
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\nwritten {out_path}")
