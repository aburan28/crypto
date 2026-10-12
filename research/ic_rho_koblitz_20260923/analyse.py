#!/usr/bin/env python3
"""Ledger §19: grade the three targets declared in §19.1 and produce every
table §19 and the scoreboard quote, from the frozen files only.

    python3 research/ic_rho_koblitz_20260923/analyse.py

Writes analysis.json beside this file and prints the tables.  Nothing here
runs a measurement: the batch runs are batch/*.json, the prices are
prices/*.json, and every index-calculus figure is read from the frozen
docs/ic/runs file the scoreboard cites for it.
"""
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(ROOT, "docs", "ic", "runs")

# Two-sided 95 % Student t quantiles by degrees of freedom.
T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306,
        9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131}
TARGET2 = (0.09, 0.35)
TARGETS = 32


def load(*parts):
    with open(os.path.join(*parts)) as f:
        return json.load(f)


def mean(v):
    return sum(v) / len(v)


def sd(v):
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0


def t_interval(v):
    h = T975[len(v) - 1] * sd(v) / math.sqrt(len(v))
    return [mean(v) - h, mean(v) + h]


def bootstrap_ratio(num, den, draws=20000, seed=19):
    """95 % percentile interval of mean(num)/mean(den), batches resampled."""
    rng = random.Random(seed)
    out = []
    for _ in range(draws):
        a = mean([rng.choice(num) for _ in num])
        b = mean([rng.choice(den) for _ in den])
        out.append(a / b)
    out.sort()
    return [out[int(0.025 * draws)], out[int(0.975 * draws) - 1]]


# ── (a) batch rho ──────────────────────────────────────────────────

batch = {}
for n in (41, 53, 61):
    doc = load(HERE, "batch", f"k0n{n}.json")
    assert doc["operation"] == "rho-batch" and len(doc["curves"]) == 1
    c = doc["curves"][0]
    sizes = {}
    verified = targets = 0
    for s in c["sizes"]:
        vals = [run["s_per_target"] for run in s["runs"]]
        per = [t for run in s["runs"] for t in run["per_target"]]
        ok = [t for t in per if t["verified"] and t["recovered"] == t["planted"]]
        targets += len(per)
        verified += len(ok)
        canon = sum(run["counters"]["canonicalisations_uncharged"] for run in s["runs"])
        walk_ops = sum(run["counters"]["walk_operations"] for run in s["runs"])
        sizes[s["k"]] = {
            "k": s["k"], "batches": len(vals), "targets": len(per), "verified": len(ok),
            "s_per_target": vals, "mean": mean(vals), "sd": sd(vals), "ci95": t_interval(vals),
            "walk_s_per_target": mean([run["s_walk_per_target"] for run in s["runs"]]),
            "batch_law": s["batch_law"],
            "solved_on_own_trail": s["solved_on_own_trail"],
            "solved_on_an_earlier_trail": s["solved_on_an_earlier_trail"],
            "canonicalisations_per_walk_operation": canon / walk_ops,
            "fruitless_cycles": sum(run["counters"]["cycles_length_2"] + run["counters"]["cycles_length_4"]
                                    + run["counters"]["cycles_other_length"] for run in s["runs"]),
            "walks_capped": sum(run["counters"]["walks_capped"] for run in s["runs"]),
        }
    one = sizes[1]
    for k, s in sizes.items():
        s["over_k1"] = s["mean"] / one["mean"]
        s["over_k1_ci95"] = bootstrap_ratio(s["s_per_target"], one["s_per_target"]) if k != 1 else [1.0, 1.0]
        s["over_batch_law"] = s["over_k1"] / s["batch_law"]
        s["over_floor"] = s["mean"] / c["floor_s_single"]
    batch[n] = {
        "instance": c["instance"], "r": c["r"], "sqrt_r": math.sqrt(c["r"]), "log2_r": c["log2_r"],
        "automorphisms": c["automorphisms"], "floor_s_single": c["floor_s_single"],
        "targets": targets, "verified": verified, "sizes": sizes,
        "binary_blake3": doc["software"]["binary_blake3"],
    }

target1 = {
    "declared": "every target of every batch is recovered and verified",
    "targets": sum(b["targets"] for b in batch.values()),
    "verified": sum(b["verified"] for b in batch.values()),
}
target1["met"] = target1["targets"] == target1["verified"]

target2 = {"declared": f"per-target cost at k = 32 over k = 1 in [{TARGET2[0]}, {TARGET2[1]}]", "curves": {}}
for n, b in batch.items():
    s = b["sizes"][32]
    target2["curves"][n] = {"ratio": s["over_k1"], "ci95": s["over_k1_ci95"], "batch_law": s["batch_law"],
                            "met": TARGET2[0] <= s["over_k1"] <= TARGET2[1]}
target2["met"] = all(c["met"] for c in target2["curves"].values())

# ── (b) and (c) prices ─────────────────────────────────────────────

prices = {}
for name in sorted(os.listdir(os.path.join(HERE, "prices"))):
    d = load(HERE, "prices", name)
    prices[name[:-5]] = d
step = {}
for n, key in ((41, "n41-F15744"), (53, "n53"), (61, "n61")):
    d = prices[key]
    step[n] = {
        "source": f"prices/{key}.json",
        "unit_ns": d["unit_ns"]["median"],
        "canonical": d["b_canonical_step_units"],
        "canonical_range": [1 + d["b_canonical_canonicalisation_units"]["min"],
                            1 + d["b_canonical_canonicalisation_units"]["max"]],
        "bailey": d["b_bailey_step_units"],
        "bailey_range": [1 + d["b_bailey_overhead_units"]["min"], 1 + d["b_bailey_overhead_units"]["max"]],
        "add_pairwise": d["diagnostic_add_pairwise_units"]["median"],
        "affine_add": d["diagnostic_affine_add_units"]["median"],
        "canonicalisation_in_affine_units": d["b_canonical_canonicalisation_units"]["median"]
        / d["diagnostic_affine_add_units"]["median"],
    }
build_price = {}
for key, d in prices.items():
    if d.get("c_build"):
        c = d["c_build"]
        build_price[(d["n"], c["base_points"])] = {
            "source": f"prices/{key}.json", "stored_pairs": c["stored_pairs"], "orbits": c["signed_orbits"],
            "units_per_stored_pair": c["units_per_stored_pair"]["median"],
            "range": [c["units_per_stored_pair"]["min"], c["units_per_stored_pair"]["max"]],
            "declared": (d["n"], c["base_points"]) == (41, 15744),
        }

# ── The figures the scoreboard quotes, from their frozen files ─────

figures = []


def figure(panel, row, n, points, tier, total, build, descent, quoted_s, quoted_vs, rho_s, source, note=None,
           lower_bound=False):
    figures.append({
        "panel": panel, "row": row, "n": n, "points": points, "tier": tier, "total": total,
        "build_counted": build, "descent": descent, "quoted_s": quoted_s, "quoted_vs_rho": quoted_vs,
        "quoted_rho_s": rho_s, "source": source, "note": note,
        # S left selection and the linear algebra unpriced: every ratio is a lower bound.
        "lower_bound": lower_bound,
    })


sp = load(RUNS, "koblitz-select-packed-20260922.json")
rho = sp["boundaries"]["rho_reference_S"]["value"]
for r in sp["ledger"]["rows"]:
    figure("koblitz-select-packed-20260922", r["variant"], 41, 15744, "folded", r["total_adds"], r["build"],
           r["descent"], r["S"], r["vs_rho"], rho, "koblitz-select-packed-20260922.json ledger.rows")

aim = load(RUNS, "koblitz-collection-aim-20260922.json")
rho = aim["boundaries"]["rho_reference_S"]["value"]
# The A/B arms are the ladder's last rows; the ledger rows carry their totals.
rows = aim["ledger"]["rows"]
# descent: recorded in select-packed for the same configuration (54 adds for
# 32 targets); the 5,248-point rows share the phase-prices rung's 864.
for r in rows:
    pts = r["points"]
    descent = 864 if pts == 5248 else (54 if pts == 15744 else None)
    figure("koblitz-collection-aim-20260922", r["variant"], 41, pts, "folded", r["total_adds"],
           None, descent, r["S"], r["vs_rho"], rho, "koblitz-collection-aim-20260922.json ledger.rows",
           None if descent is not None else "descent not recorded; cold figure includes every target's descent")

pp = load(RUNS, "koblitz-phase-prices-20260921.json")
for r in pp["part_1_phase_prices"]["rungs"]:
    if not r.get("priced"):
        continue
    ph = r["phases_group_additions"]
    figure("koblitz-phase-prices-20260921", r["rung"], r["degree"], r["points"], r["tier"],
           r["total_group_additions"], ph["pair_table_build"], ph["descent_probes"], r["S"], r["over_rho"],
           r["rho_S"], "koblitz-phase-prices-20260921.json part_1_phase_prices.rungs")

cw = load(RUNS, "koblitz-collection-window-20260921.json")
# The page shows the re-priced values; its n = 41 rows and the holdout are
# quoted against 0.1665 and 0.2172 (see the page), so the quoted ratio is
# recomputed from the page's reference, not the file's.
page_rho = {41: 0.1665, 53: 0.2172}
for v in cw["variants"]:
    ph = v["phases_group_additions_repriced_20260922"]
    s = v["S_repriced_20260922"]
    figure("koblitz-collection-window-20260921", v["variant"], v["degree"], v["points"], v["tier"],
           v["total_group_additions_repriced_20260922"], ph["pair_table_build"], ph["descent_probes"], s,
           round(s / page_rho[v["degree"]], 2), page_rho[v["degree"]],
           "koblitz-collection-window-20260921.json variants (*_repriced_20260922)",
           "file quotes rho_S = %s here" % v["rho_S"] if v["rho_S"] != page_rho[v["degree"]] else None)

pv = load(RUNS, "koblitz-probe-volume-20260921.json")
for r in pv["rungs"]:
    if not r.get("tiers"):
        continue
    for tier in ("compact", "folded"):
        t = r["tiers"][tier]
        figure("koblitz-probe-volume-20260921", f"{r['rung']} {tier}", r["degree"], r["points"], tier,
               t["total_adds"], t["build_adds"], t["descent_adds"], t["S"], t["over_rho"], r["rho_S"],
               "koblitz-probe-volume-20260921.json rungs[].tiers",
               "S is a lower bound: selection and the linear algebra unpriced", lower_bound=True)

tc = load(RUNS, "koblitz-tier-crossover-20260921.json")
for w in tc["widths"]:
    for tier in ("full", "compact", "folded"):
        t = w["tiers"][tier]
        figure("koblitz-tier-crossover-20260921", f"{w['points']} points {tier}", 61, w["points"], tier,
               t["total_adds"], t["build_adds"], t["descent_adds"], t["S_lower_bound"], t["over_rho"],
               tc["boundaries"]["rho_S"], "koblitz-tier-crossover-20260921.json widths[].tiers",
               "S is a lower bound: selection and the linear algebra unpriced", lower_bound=True)

# ── The re-read ────────────────────────────────────────────────────


def stored_pairs(n, points):
    bp = build_price.get((n, points))
    return bp["stored_pairs"] if bp else None


for f in figures:
    b = batch[f["n"]]
    sqrt_r = b["sqrt_r"]
    s = f["total"] / (TARGETS * sqrt_r)
    assert abs(s - f["quoted_s"]) <= 0.0006 * max(1.0, f["quoted_s"]), (f["row"], s, f["quoted_s"])
    b32, b1 = b["sizes"][32]["mean"], b["sizes"][1]["mean"]
    st = step[f["n"]]
    f["s"] = s
    f["reread"] = s / b32
    f["with_b_canonical"] = s / (b32 * st["canonical"])
    f["with_b_bailey_model"] = s / (b32 * st["bailey"])
    bp = build_price.get((f["n"], f["points"])) if f["tier"] == "folded" else None
    if bp:
        # The build's measured cost replaces its count: every stored pair
        # at the measured price, instead of one addition per counted pair.
        build_counted = f["build_counted"] if f["build_counted"] is not None else bp["stored_pairs"]
        total_c = f["total"] - build_counted + bp["stored_pairs"] * bp["units_per_stored_pair"]
        f["build_price_source"] = bp["source"]
        f["s_with_c"] = total_c / (TARGETS * sqrt_r)
        f["with_c"] = f["s_with_c"] / b32
        f["with_b_and_c"] = f["s_with_c"] / (b32 * st["canonical"])
    else:
        total_c = None
        f["build_price_source"] = None
        f["s_with_c"] = f["with_c"] = f["with_b_and_c"] = None
    descent = f["descent"] or 0
    cold_total = f["total"] - descent * (TARGETS - 1) / TARGETS
    f["cold_s"] = cold_total / sqrt_r
    f["cold"] = f["cold_s"] / b1
    f["cold_with_b_and_c"] = (
        (total_c - descent * (TARGETS - 1) / TARGETS) / sqrt_r / (b1 * st["canonical"]) if total_c else None
    )

target3 = {
    "declared": "every 32-target Koblitz figure on the page re-read against batch rho at k = 32 on its own curve, "
                "(b) and (c) beside it, the one-target figure beside those",
    "figures": len(figures),
    "curves_measured": sorted(batch),
    "met": all(f["n"] in batch for f in figures),
}

# ── Output ─────────────────────────────────────────────────────────


def fx(v, d=2, bound=False):
    if v is None:
        return "—"
    g = "≥ " if bound else ""
    if v >= 100:
        return f"{g}{v:,.0f}×"
    if v >= 10:
        return f"{g}{v:.1f}×"
    return f"{g}{v:.{d}f}×"


lines = ["| curve | k | batches | S per target (95 % CI) | over k = 1 (95 % CI) | batch law | over floor | own / earlier trail | ok |",
         "|:--|--:|--:|--:|--:|--:|--:|--:|:--|"]
for n, b in batch.items():
    for k, s in sorted(b["sizes"].items()):
        lines.append(
            f"| {b['instance']} | {k} | {s['batches']} | {s['mean']:.4f} [{s['ci95'][0]:.4f}, {s['ci95'][1]:.4f}] | "
            f"{s['over_k1']:.3f} [{s['over_k1_ci95'][0]:.3f}, {s['over_k1_ci95'][1]:.3f}] | {s['batch_law']:.3f} | "
            f"{s['over_floor']:.3f} | {s['solved_on_own_trail']} / {s['solved_on_an_earlier_trail']} | "
            f"{s['verified']}/{s['targets']} |")
batch_md = "\n".join(lines)

lines = ["| curve | unit (ns) | canonical step | Bailey step | add_pairwise | affine add | canonicalisation in affine units |",
         "|:--|--:|--:|--:|--:|--:|--:|"]
for n, st in step.items():
    lines.append(f"| n = {n} | {st['unit_ns']:.1f} | {st['canonical']:.2f} [{st['canonical_range'][0]:.2f}, {st['canonical_range'][1]:.2f}] | "
                 f"{st['bailey']:.2f} [{st['bailey_range'][0]:.2f}, {st['bailey_range'][1]:.2f}] | {st['add_pairwise']:.2f} | "
                 f"{st['affine_add']:.1f} | {st['canonicalisation_in_affine_units']:.3f} |")
step_md = "\n".join(lines)

lines = ["| base | orbits | stored pairs | units per stored pair | range | declared |", "|:--|--:|--:|--:|--:|:--|"]
for (n, pts), bp in sorted(build_price.items()):
    lines.append(f"| n = {n}, F = {pts:,} | {bp['orbits']} | {bp['stored_pairs']:,} | {bp['units_per_stored_pair']:.2f} | "
                 f"[{bp['range'][0]:.2f}, {bp['range'][1]:.2f}] | {'yes' if bp['declared'] else 'no'} |")
build_md = "\n".join(lines)

lines = ["| panel | row | n | quoted | re-read (k = 32) | + (b) canonical | + (b) Bailey, model | + (c) | + (b) and (c) | cold | cold + (b), (c) |",
         "|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
for f in figures:
    lb = f["lower_bound"]
    lines.append(f"| {f['panel'].replace('koblitz-', '')} | {f['row']} | {f['n']} | {fx(f['quoted_vs_rho'], bound=lb)} | "
                 f"**{fx(f['reread'], bound=lb)}** | {fx(f['with_b_canonical'], bound=lb)} | "
                 f"{fx(f['with_b_bailey_model'], bound=lb)} | {fx(f['with_c'], bound=lb)} | "
                 f"{fx(f['with_b_and_c'], bound=lb)} | {fx(f['cold'], 1, lb)} | {fx(f['cold_with_b_and_c'], 1, lb)} |")
reread_md = "\n".join(lines)

analysis = {
    "what_this_is": "Ledger section 19: the Koblitz collection thread's 32-target figures re-read against batch rho "
                    "at k = 32 on the same curve, with the measured step and build prices beside them.",
    "targets": {"1": target1, "2": target2, "3": target3},
    "batch": {str(n): b for n, b in batch.items()},
    "step_price": {str(n): s for n, s in step.items()},
    "build_price": {f"n{n}-F{p}": bp for (n, p), bp in build_price.items()},
    "figures": figures,
    "markdown": {"batch": batch_md, "step": step_md, "build": build_md, "reread": reread_md},
}
with open(os.path.join(HERE, "analysis.json"), "w") as f:
    json.dump(analysis, f, indent=1, sort_keys=True)
    f.write("\n")

print(batch_md, "\n")
print(step_md, "\n")
print(build_md, "\n")
print(reread_md, "\n")
for k, t in analysis["targets"].items():
    print(f"target {k}: {'MET' if t['met'] else 'NOT MET'} — {t['declared']}")
print(json.dumps(target2["curves"], indent=1))
