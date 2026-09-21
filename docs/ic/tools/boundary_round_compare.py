#!/usr/bin/env python3
"""Pair the before and after rows of a boundary-ledger run and save the
baseline/candidate comparison AGENTS.md §8 asks for.

usage: boundary_round_compare.py <round.json> [--baseline <round1.json>] [--out <comparison.json>] [--holdout <holdout.json>]

Every Round-2 row is a rung of a cumulative ledger named by suffix:
`<first-round row>` → `_negfold` → `_frobfold` → `_walk`, with `_balanced`
for the resized Koblitz base and `mitm_m2_…` rows for the two-summand count.
For each candidate row the script finds its previous rung (the same name
with the last suffix stripped, or the first-round best of the regime when
that rung was never run) and reports

    speedup = baseline_total_operations / candidate_total_operations

per repeat (the runs are paired: same curve, same target, same seed) and
as the mean, together with S, the ratios to rho and to the floor, and the
correctness column.  It also checks that the first-round rows of the new
run reproduce the frozen Round-1 counts exactly (same seeds, same code
path), which is what makes the pairing a comparison rather than two runs.
"""
import json, math, sys
from collections import OrderedDict

def arg(flag, default=None):
    if flag in sys.argv:
        return sys.argv[sys.argv.index(flag) + 1]
    return default

run_path = sys.argv[1]
run = json.load(open(run_path))
baseline = json.load(open(arg("--baseline"))) if arg("--baseline") else None
holdout = json.load(open(arg("--holdout"))) if arg("--holdout") else None
out_path = arg("--out")

SUFFIXES = ("_walk", "_frobfold", "_negfold", "_balanced")
ROUND2 = ("_negfold", "_frobfold", "_walk", "_balanced")
CLASS = {
    "_negfold": "engineering",
    "_frobfold": "engineering",
    "_walk": "engineering",
    "_balanced": "engineering",
}

def is_round2(name):
    return any(s in name for s in ROUND2)

def previous_rung(name):
    """The row this one was built on: strip the last Round-2 suffix."""
    for s in SUFFIXES:
        if name.endswith(s):
            return name[: -len(s)]
    return None

def rows_of(inst):
    seen = OrderedDict()
    for v in inst["variants"]:
        seen.setdefault(v["name"], []).append(v)
    return seen

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else float("nan")

def f(v, d=3):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "—"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    if abs(v) >= 100:
        return f"{v:.0f}"
    if abs(v) >= 10:
        return f"{v:.1f}"
    return f"{v:.{d}f}"

def compare_instance(inst):
    rows = rows_of(inst)
    first_round = {n: rs for n, rs in rows.items() if not is_round2(n)}
    if not first_round:
        return []
    best_before_name = min(first_round, key=lambda n: mean([r["s"] for r in first_round[n]]))
    out = []
    for name, runs in rows.items():
        if not is_round2(name):
            continue
        prev = previous_rung(name)
        note = ""
        if prev not in rows:
            # A rung that was never run (e.g. the two-summand row's
            # non-walk version at large n, or a balanced row whose narrow
            # counterpart differs in m): compare with the first-round best.
            prev = best_before_name
            note = "previous rung not run; against the first-round best"
        base_runs = rows[prev]
        if runs[0]["summands"] != base_runs[0]["summands"]:
            note = (note + "; " if note else "") + "summand count changed"
        if runs[0]["signed_points"] != base_runs[0]["signed_points"]:
            note = (note + "; " if note else "") + "base changed"
        k = min(len(runs), len(base_runs))
        paired = [base_runs[i]["total_gae"] / runs[i]["total_gae"] for i in range(k)]
        best_runs = rows[best_before_name]
        vs_best = [best_runs[i]["total_gae"] / runs[i]["total_gae"] for i in range(min(len(runs), len(best_runs)))]
        lever = next((s for s in SUFFIXES if name.endswith(s)), "")
        out.append(OrderedDict([
            ("regime", inst["regime"]),
            ("instance", inst["curve"]["name"]),
            ("log2_r", round(inst["log2_r"], 2)),
            ("candidate", name),
            ("baseline", prev),
            ("lever", lever.strip("_")),
            ("class", CLASS.get(lever, "engineering")),
            ("note", note),
            ("m", runs[0]["summands"]),
            ("signed_points", runs[0]["signed_points"]),
            ("columns", runs[0]["columns"]),
            ("repeats_paired", k),
            ("speedup_per_repeat", [round(x, 4) for x in paired]),
            ("speedup_mean", round(mean(paired), 4)),
            ("speedup_min", round(min(paired), 4)),
            ("speedup_max", round(max(paired), 4)),
            ("speedup_vs_first_round_best", round(mean(vs_best), 4)),
            ("first_round_best", best_before_name),
            ("S_baseline", round(mean([r["s"] for r in base_runs]), 4)),
            ("S_candidate", round(mean([r["s"] for r in runs]), 4)),
            ("S_over_rho_baseline", round(mean([r["s"] for r in base_runs]) / inst["rho_s_mean"], 4)),
            ("S_over_rho_candidate", round(mean([r["s"] for r in runs]) / inst["rho_s_mean"], 4)),
            ("S_over_floor_candidate", round(mean([r["s"] for r in runs]) / inst["floor_s"], 4)),
            ("yield_over_ceiling_candidate", round(mean([r["yield_over_ceiling"] for r in runs]), 4)),
            ("yield_over_ceiling_exact_candidate", round(mean([r["yield_over_ceiling_exact"] for r in runs]), 4)),
            ("trials_candidate", round(mean([r["trials"] for r in runs]), 1)),
            ("verified_candidate", all(r["verified"] for r in runs)),
            ("verified_baseline", all(r["verified"] for r in base_runs)),
            ("rho_verified", inst["rho_verified_all"]),
        ]))
    return out

def reproduction_check(new, old):
    """First-round rows of the new run against the frozen Round-1 file:
    the native counts must be identical (same seeds, same code path)."""
    checks = []
    old_by = {(i["regime"], i["curve"]["name"]): i for i in old["ledger"]["instances"]}
    for inst in new["ledger"]["instances"]:
        key = (inst["regime"], inst["curve"]["name"])
        if key not in old_by:
            continue
        old_rows = rows_of(old_by[key])
        for name, runs in rows_of(inst).items():
            if is_round2(name) or name not in old_rows:
                continue
            o = old_rows[name]
            same = len(o) == len(runs) and all(
                a["trials"] == b["trials"]
                and a["relations_found"] == b["relations_found"]
                and a["relations"]["group_ops"] == b["relations"]["group_ops"]
                and a["factor_base"]["group_ops"] == b["factor_base"]["group_ops"]
                and a["linear_algebra"]["native"].get("row_ops") == b["linear_algebra"]["native"].get("row_ops")
                and a["recovered"] == b["recovered"]
                for a, b in zip(o, runs)
            )
            checks.append(OrderedDict([("regime", inst["regime"]), ("instance", inst["curve"]["name"]), ("variant", name), ("counts_identical", same)]))
    return checks

comparison = []
for inst in run["ledger"]["instances"]:
    comparison.extend(compare_instance(inst))

repro = reproduction_check(run, baseline) if baseline else []
holdout_rows = []
if holdout:
    for inst in holdout["ledger"]["instances"]:
        holdout_rows.extend(compare_instance(inst))

# ── Markdown ──
print("| regime | instance | log2 r | candidate | baseline | lever | class | m | |F| | K | speedup per repeat | speedup mean | vs first-round best | S before → after | vs rho after | vs floor after | y/c exact | ok |")
print("|:--|:--|--:|:--|:--|:--|:--|--:|--:|--:|:--|--:|--:|:--|--:|--:|--:|:--|")
for c in comparison:
    print(f"| {c['regime']} | {c['instance']} | {c['log2_r']:.1f} | {c['candidate']} | {c['baseline']} | {c['lever']} | {c['class']} | {c['m']} | {c['signed_points']} | {c['columns']} | {', '.join(f(x) for x in c['speedup_per_repeat'])} | {f(c['speedup_mean'])} | {f(c['speedup_vs_first_round_best'])} | {f(c['S_baseline'])} → {f(c['S_candidate'])} | {f(c['S_over_rho_candidate'])}× | {f(c['S_over_floor_candidate'])}× | {f(c['yield_over_ceiling_exact_candidate'])} | {'✓' if c['verified_candidate'] and c['verified_baseline'] and c['rho_verified'] else '✗'} |")
if repro:
    bad = [r for r in repro if not r["counts_identical"]]
    print(f"\nFirst-round rows reproduced against the frozen Round-1 file: {len(repro) - len(bad)} of {len(repro)} identical in every native count.")
    for r in bad:
        print(f"  differs: {r}")
if holdout_rows:
    print("\nHoldout (fresh seed):\n")
    print("| regime | instance | candidate | baseline | speedup per repeat | speedup mean | S after | vs rho after | ok |")
    print("|:--|:--|:--|:--|:--|--:|--:|--:|:--|")
    for c in holdout_rows:
        print(f"| {c['regime']} | {c['instance']} | {c['candidate']} | {c['baseline']} | {', '.join(f(x) for x in c['speedup_per_repeat'])} | {f(c['speedup_mean'])} | {f(c['S_candidate'])} | {f(c['S_over_rho_candidate'])}× | {'✓' if c['verified_candidate'] and c['verified_baseline'] and c['rho_verified'] else '✗'} |")

if out_path:
    doc = OrderedDict([
        ("what", "Baseline/candidate comparison of the Round-2 boundary-ledger run: every Round-2 row paired with the rung it was built on, on identical curves, subgroups, factor bases, targets and seeds; speedup = baseline_total_operations / candidate_total_operations in group-addition equivalents with every phase charged"),
        ("unit", run["ledger"]["unit"]),
        ("run", run_path),
        ("baseline_round1", arg("--baseline")),
        ("holdout", arg("--holdout")),
        ("regression_suite_note", "The frozen WDSat solver regression suite (research/index_calculus_baseline_20260914/regression/) measures a SAT-solver stage on 60 fixed inputs and is not an adapter for other pipelines; the levers here (pair-table folds, target generation, exact ceiling, base sizing) touch no solver stage, so the matched full-DLP comparison in this file is the §8 evidence, under the parent accounting contract's exclusive per-phase charging."),
        ("all_candidates_verified", all(c["verified_candidate"] for c in comparison)),
        ("all_baselines_verified", all(c["verified_baseline"] for c in comparison)),
        ("first_round_rows_reproduced", all(r["counts_identical"] for r in repro) if repro else None),
        ("reproduction_checks", repro),
        ("comparison", comparison),
        ("holdout_comparison", holdout_rows),
    ])
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\nwritten {out_path}")
