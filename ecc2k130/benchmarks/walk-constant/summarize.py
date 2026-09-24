"""Print WALK-CONSTANT.md's tables from the frozen run files in this directory.

Every number in the note's tables comes from here: matrix-v2.jsonl (the
emulation, examples/ecc2k130_walk_constant.rs, table and adding walks on the
device's cycle rule), matrix-v1.jsonl (its sigma rows; its table and adding
rows used a previous-class rule only and are superseded), device-v2.jsonl and
device-v1.jsonl (ecc2k130's own reference walks, src/walkconstant.cpp).
The sigma walk at n = 19 is degenerate and excluded from the aggregates: two of
its eight multipliers fall in one class (1 + s^10 = s^-9 (1 + s^9) when
n = 19; every n >= 21 has eight distinct ones), see sigma_classes.py.  The
cost-to-solve table combines the aggregates with the paired rates of
ITERATION-FUNCTION.md section 6 and ONE-BLOCK-GEOMETRY.md section 1 and the loss
model of trap_cost.py.
Run: python3 summarize.py
"""
import json
from math import log2, sqrt

import trap_cost


def load(name):
    try:
        with open(name) as f:
            return [dict(json.loads(line), source=name) for line in f if line.strip()]
    except FileNotFoundError:
        return []


v1 = load("matrix-v1.jsonl")
v2 = load("matrix-v2.jsonl")
d1 = load("device-v1.jsonl")
d2 = load("device-v2.jsonl")

rows = [r for r in v1 if r["walk"] == "sigma"] + v2 + d1 + d2


def label(r):
    walk = r["walk"].replace("device-", "")
    where = "device" if r["walk"].startswith("device-") else "emulation"
    fixed = ", one mapping" if r.get("fixed_mapping") else ""
    return f"{walk} | {r['branches']} | {r['dist']}{fixed} | {where}"


def model(r):
    n, s2 = r["n"], r["sum_p2"]
    injective = 1 / sqrt(1 - s2)
    class_frame = 1 / sqrt(1 - s2 / (2 * n))
    return class_frame if "sigma" not in r["walk"] else injective


print("| n | walk | H | branches | harness | W | trials | c | ± | model | c / model | source |")
print("|---:|---|---:|---|---|---:|---:|---:|---:|---:|---:|---|")
for r in sorted(rows, key=lambda r: (r["n"], "sigma" in r["walk"], r["walk"], r["branches"], r["dist"],
                                     r.get("fixed_mapping", False), r["walks"], r["seed"])):
    if r.get("fixed_mapping"):
        continue
    m = model(r)
    print(f"| {r['n']} | {label(r)} | {r['walks']} | {r['trials']:,} | {r['c']:.4f} | {r['c_se']:.4f} | "
          f"{m:.4f} | {r['c'] / m:.4f} | {r['source']} |")

print()
print("One mapping at a time (W = 16):")
print("| n | walk | mappings (seed: c) | mean | spread (sd) | per-mapping ± |")
print("|---:|---|---|---:|---:|---:|")
spread = {}
for key in sorted({(r["n"], r["walk"]) for r in rows if r.get("fixed_mapping")}):
    sel = [r for r in rows if r.get("fixed_mapping") and (r["n"], r["walk"]) == key]
    cs = [r["c"] for r in sel]
    m = sum(cs) / len(cs)
    sd = sqrt(sum((c - m) ** 2 for c in cs) / (len(cs) - 1)) if len(cs) > 1 else float("nan")
    se = sum(r["c_se"] for r in sel) / len(sel)
    spread[key] = (m, sd, se, len(cs))
    listing = ", ".join(f"{r['seed']}: {r['c']:.4f}" for r in sorted(sel, key=lambda r: r["seed"]))
    print(f"| {key[0]} | {key[1]} ecc2k130 H = 8 | {listing} | {m:.4f} | {sd:.4f} | {se:.4f} |")

print()
print("Fruitless returns, measured against the leading-order prediction:")
print("| n | walk | H | branches | harness | steps | pairwise seen | expected | τ-relation seen | expected |")
print("|---:|---|---:|---|---|---:|---:|---:|---:|---:|")
for r in sorted(v2 + d2, key=lambda r: (r["n"], r["walk"], r["branches"], r["dist"])):
    if "sigma" in r["walk"] or r.get("fixed_mapping") or "relation" not in r:
        continue
    steps = r["steps"]
    pair = sum(r["fruitless"].values())
    rel = sum(r["relation"].values())
    print(f"| {r['n']} | {label(r)} | {steps:,} | {pair} | {r['fruitless_predicted'] * steps:.1f} | "
          f"{rel} | {r['relation_predicted'] * steps:.1f} |")
tot = {}
for r in v2 + d2:
    if "sigma" in r["walk"] or "relation" not in r:
        continue
    for k, seen, pred in (("pairwise", sum(r["fruitless"].values()), r["fruitless_predicted"] * r["steps"]),
                          ("relation", sum(r["relation"].values()), r["relation_predicted"] * r["steps"])):
        a, b = tot.get(k, (0, 0.0))
        tot[k] = (a + seen, b + pred)
for k, (seen, pred) in tot.items():
    print(f"all rows, {k}: {seen} seen, {pred:.1f} expected, ratio {seen / pred:.3f} ± {sqrt(seen) / pred:.3f}")
for k in ("pairwise", "relation"):
    seen = pred = 0
    for r in v2 + d2:
        if "sigma" in r["walk"] or "relation" not in r or r["n"] < 37:
            continue
        key = "fruitless" if k == "pairwise" else "relation"
        seen += sum(r[key].values())
        pred += r[key + "_predicted"] * r["steps"]
    if pred:
        print(f"rows at n >= 37, {k}: {seen} seen, {pred:.1f} expected, ratio {seen / pred:.3f} ± {sqrt(seen) / pred:.3f}")


def pooled(pick):
    sel = [r for r in rows if pick(r) and not r.get("fixed_mapping")]
    w = [1 / r["c_se"] ** 2 for r in sel]
    c = sum(wi * r["c"] for wi, r in zip(w, sel)) / sum(w)
    return c, sqrt(1 / sum(w)), sel


def hashed131(walk):
    return lambda r: (r["walk"].replace("device-", "") == walk and r["dist"] == "ecc2k130"
                      and r["branches"] == 8 and not (walk == "sigma" and r["n"] < 21))


c_table, se_table, sel_t = pooled(hashed131("table"))
print()
print(f"table walk, ECC2K-130 branch distribution, H = 8, pooled: c = {c_table:.4f} ± {se_table:.4f} "
      f"over {len(sel_t)} rows (n = {sorted(set(r['n'] for r in sel_t))})")
# The sigma constant falls with the degree, so it is not pooled: at n = 131
# it is bracketed by its first-order value (the limit, from above) and the
# largest degree measured.  That bracket is an extrapolation.
by_n = {}
for r in rows:
    if hashed131("sigma")(r) and not r.get("fixed_mapping") and r["walk"] == "sigma":
        by_n.setdefault(r["n"], []).append(r)
trend = []
for n in sorted(by_n):
    w = [1 / r["c_se"] ** 2 for r in by_n[n]]
    c = sum(wi * r["c"] for wi, r in zip(w, by_n[n])) / sum(w)
    trend.append((n, c, sqrt(1 / sum(w))))
print("sigma walk, ECC2K-130 branch distribution, emulation, by degree: " +
      ", ".join(f"n = {n}: {c:.4f} ± {se:.4f}" for n, c, se in trend))
first_order = 1 / sqrt(1 - sum(x * x for x in trap_cost.branch_probabilities(131, 8)))
sigma_lo, sigma_hi = first_order, trend[-1][1]
ratio_lo, ratio_hi = c_table / sigma_hi, c_table / sigma_lo
print(f"sigma at n = 131, bracketed (extrapolation): [{sigma_lo:.4f}, {sigma_hi:.4f}]; "
      f"table / sigma in iterations: {ratio_lo:.4f} - {ratio_hi:.4f}")
print(f"sigma expected work on completed trails: 2^60.809 x [{sigma_lo:.4f}, {sigma_hi:.4f}] = "
      f"2^[{60.809 + log2(sigma_lo):.3f}, {60.809 + log2(sigma_hi):.3f}]; with the current guard (x1.185 at dpWeight 32): "
      f"2^[{60.809 + log2(sigma_lo * 1.185):.3f}, {60.809 + log2(sigma_hi * 1.185):.3f}]")

# Cost to solve, table over sigma, at the live dpWeight = 32.
PAIRED = (("automatic workers", 14.41, 16.56), ("audited 385k", 14.98, 16.35), ("one-block session 256x2", 14.00, 15.84))
th = trap_cost.theta(131, 32)
p8 = trap_cost.branch_probabilities(131, 8)
s2, s4 = sum(x * x for x in p8), sum(x ** 4 for x in p8)
pair, rel = 4 * (s2 / 262) ** 3, 24 * s4 / 262 ** 3
loss_sigma_2_30 = trap_cost.overhead(0.0, th, 2 ** 30)
# The measured fruitless rates run below the leading-order count; price
# with the count and with the count scaled to the measurements.
seen_all = sum(sum(r["fruitless"].values()) + sum(r["relation"].values())
               for r in v2 + d2 if r["walk"] in ("table", "device-table"))
pred_all = sum((r["fruitless_predicted"] + r["relation_predicted"]) * r["steps"]
               for r in v2 + d2 if r["walk"] in ("table", "device-table"))
scale = seen_all / pred_all
print()
print(f"fruitless returns on all table rows: {seen_all} seen, {pred_all:.0f} predicted, scale {scale:.3f}")
print("cost to solve, table / sigma, dpWeight = 32, H = 8 (the count; the count x measured scale):")
for label, table_loss, sigma_loss in (
    ("as built, both at maxIters = 2^30",
     lambda k: trap_cost.overhead(k * (pair + rel), th, 2 ** 30), loss_sigma_2_30),
    ("as built, each at its best guard", lambda k: trap_cost.best_guard(k * (pair + rel), th)[0], 1.0),
    ("rule also refusing tau-relations, best guards", lambda k: trap_cost.best_guard(k * pair, th)[0], 1.0),
    ("every fruitless cycle caught and escaped, best guards", lambda k: 1.0, 1.0),
):
    losses = [table_loss(k) for k in (1.0, scale)]
    vals = [q * (rs / rt) * lt / sigma_loss for _, rs, rt in PAIRED for q in (ratio_lo, ratio_hi) for lt in losses]
    print(f"  {label}: table loss x{losses[0]:.3f} (x{losses[1]:.3f} scaled), sigma loss x{sigma_loss:.3f}; "
          f"table / sigma = {min(vals):.2f} - {max(vals):.2f}")


# ---------------------------------------------------------------------------
# Round 2 (WALK-CONSTANT.md section 11): the table walk under the extended
# cycle rule, v2.  matrix-v3 is the emulation under v2, device-v3 the device's
# own walk (whose rule is v2 since that section), merge.jsonl the merge
# parting under both rules.  Everything above is rule v1 and stays as it was.
v3 = load("matrix-v3.jsonl")
d3 = load("device-v3.jsonl")
merges = load("merge.jsonl")

print()
print("Round 2: the table walk under the extended rule (v2), against the same rows under v1:")
print("| n | branches | H | harness | W | trials | c (v2) | ± | c (v1) | ± | difference / SE | "
      "pairwise returns (v1 count) | τ-relation returns (v1 count) | rule fired per step |")
print("|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
old_rows = [r for r in v2 + d2 if r["walk"] in ("table", "device-table") and not r.get("fixed_mapping")]
for r in sorted(v3 + d3, key=lambda r: (r["n"], r["dist"], r["branches"], r["walk"])):
    if r["walk"] not in ("table", "device-table") or r.get("fixed_mapping"):
        continue
    same = [o for o in old_rows if (o["n"], o["dist"], o["branches"], o["walk"], o["walks"]) ==
            (r["n"], r["dist"], r["branches"], r["walk"], r["walks"])]
    old = same[0] if same else None
    where = "device" if r["walk"].startswith("device") else "emulation"
    pair = sum(r["fruitless"].values())
    rel = sum(r["relation"].values())
    diff = (f"{(r['c'] - old['c']) / sqrt(r['c_se'] ** 2 + old['c_se'] ** 2):+.1f}" if old else "–")
    oldc = f"{old['c']:.4f} | {old['c_se']:.4f}" if old else "– | –"
    print(f"| {r['n']} | {r['dist']} | {r['branches']} | {where} | {r['walks']} | {r['trials']:,} | {r['c']:.4f} | "
          f"{r['c_se']:.4f} | {oldc} | {diff} | {pair} ({r['fruitless_predicted'] * r['steps']:.1f}) | "
          f"{rel} ({r['relation_predicted'] * r['steps']:.1f}) | {r['cycle_rule'] / r['steps']:.2e} |")
seen_v3 = sum(sum(r["fruitless"].values()) + sum(r["relation"].values()) for r in v3 + d3
              if r["walk"] in ("table", "device-table"))
count_v3 = sum((r["fruitless_predicted"] + r["relation_predicted"]) * r["steps"] for r in v3 + d3
               if r["walk"] in ("table", "device-table"))
steps_v3 = sum(r["steps"] for r in v3 + d3 if r["walk"] in ("table", "device-table"))
print(f"formal returns under v2: {seen_v3} in {steps_v3:,} table steps, where rule v1's count predicts "
      f"{count_v3:.0f}")
c_table_v2, se_table_v2, sel_v2 = 0.0, 0.0, []
sel_v2 = [r for r in v3 + d3 if r["walk"] in ("table", "device-table") and r["dist"] == "ecc2k130"
          and r["branches"] == 8 and not r.get("fixed_mapping")]
if sel_v2:
    w = [1 / r["c_se"] ** 2 for r in sel_v2]
    c_table_v2 = sum(wi * r["c"] for wi, r in zip(w, sel_v2)) / sum(w)
    se_table_v2 = sqrt(1 / sum(w))
    print(f"table walk under v2, ECC2K-130 branches, H = 8, pooled: c = {c_table_v2:.4f} ± {se_table_v2:.4f} "
          f"over {len(sel_v2)} rows (n = {sorted(set(r['n'] for r in sel_v2))})")

print()
print("Merge parting: two walks meet at one point with different pasts; parted within 16 steps:")
print("| n | branches | rule | harness | merges | parted | rate | ± | predicted | by step (1, 2, 3, 4) |")
print("|---:|---|---|---|---:|---:|---:|---:|---:|---|")
for r in sorted(merges, key=lambda r: (r["n"], r["dist"], r["rule"], r["walk"])):
    where = "device" if r["walk"].startswith("device") else "emulation"
    steps4 = ", ".join(str(x) for x in r["parted_by_step"][:4])
    print(f"| {r['n']} | {r['dist']} | {r['rule']} | {where} | {r['merges']:,} | {r['parted']:,} | "
          f"{r['parted_rate']:.3e} | {r['parted_se']:.1e} | {r['predicted']:.3e} | {steps4} |")
ratios = {}
for r in merges:
    ratios.setdefault(r["rule"], []).append(r["parted_rate"] / r["predicted"])
for rule, xs in sorted(ratios.items()):
    print(f"rule {rule}: measured / predicted over {len(xs)} rows: {min(xs):.3f} - {max(xs):.3f}")
q131 = s2 / 262
delta_v1, delta_v2 = 2 * q131, 20 * q131
print(f"at n = 131, H = 8 (q = sum p^2 / 2n = {q131:.3e}): merges parted {delta_v1:.2%} under v1, {delta_v2:.2%} "
      f"under v2; iteration factor 1 + delta / 2 = {1 + delta_v1 / 2:.4f} and {1 + delta_v2 / 2:.4f}")

# Cost to solve with the extended rule, dpWeight 32, maxIters 2^32 (both walks
# at the new guard).  The table walk's trap loss under v2 is bounded from the
# enumeration: nothing survives at four determined tags or fewer, so the
# residual starts at order (2n)^-5; even a thousand such patterns would trap
# a fraction below 1e-5 of trails, a loss under x1.0002 at twelve trail
# lengths.  The GPU rate of the v2 kernel is not measured here; the rows use
# the v1 kernel's paired rates, so they are projections until it is.
guard = 2 ** 32
loss_sigma = trap_cost.overhead(0.0, th, guard)
trap_bound = 1000 * q131 ** 5 / th
loss_table_v2 = trap_cost.overhead(1000 * q131 ** 5, th, guard)
if c_table_v2:
    print()
    print(f"cost to solve, table (v2) / sigma, dpWeight 32, maxIters 2^32: sigma loss x{loss_sigma:.5f}, table trap "
          f"loss <= x{loss_table_v2:.5f} (bound; trapped fraction <= {trap_bound:.1e}), merge factor x{1 + delta_v2 / 2:.4f}")
    vals = [c_table_v2 * (1 + delta_v2 / 2) / cs * (rs / rt) * loss_table_v2 / loss_sigma
            for _, rs, rt in PAIRED for cs in (sigma_lo, sigma_hi)]
    print(f"  table / sigma = {min(vals):.2f} - {max(vals):.2f} (projection: the v1 kernel's paired rates)")
