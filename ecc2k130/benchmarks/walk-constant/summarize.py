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
print("One mapping at a time (n = 37, W = 16):")
print("| walk | seed | c | ± |")
print("|---|---:|---:|---:|")
for r in rows:
    if r.get("fixed_mapping"):
        print(f"| {r['walk']} {r['dist']} H = {r['branches']} | {r['seed']} | {r['c']:.4f} | {r['c_se']:.4f} |")

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
configs = (
    ("as built, both at maxIters = 2^30", trap_cost.overhead(pair + rel, th, 2 ** 30), loss_sigma_2_30),
    ("as built, each at its best guard", trap_cost.best_guard(pair + rel, th)[0], 1.0),
    ("rule also refusing tau-relations, best guards", trap_cost.best_guard(pair, th)[0], 1.0),
    ("every fruitless cycle caught and escaped, best guards", 1.0, 1.0),
)
print()
print("cost to solve, table / sigma, dpWeight = 32, H = 8:")
for label, lt, ls in configs:
    vals = [q * (rs / rt) * lt / ls for _, rs, rt in PAIRED for q in (ratio_lo, ratio_hi)]
    print(f"  {label}: table loss x{lt:.3f}, sigma loss x{ls:.3f}; table / sigma = {min(vals):.2f} - {max(vals):.2f}")
