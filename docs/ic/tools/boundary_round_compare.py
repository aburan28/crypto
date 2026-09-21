#!/usr/bin/env python3
"""Pair the before and after rows of a boundary-ledger run and save the
baseline/candidate comparison AGENTS.md §8 asks for.

usage: boundary_round_compare.py <round.json> [--baseline <round1.json>]
       [--across <previous-round.json>] [--out <comparison.json>] [--holdout <holdout.json>]

There are two pairings here and they answer different questions.

`--across` is the one AGENTS.md §8 asks for when a round changes a row
*in place* rather than adding a new one: it pairs each row of the
candidate run with the **same-named row of the previous round's run**,
on the same curve, the same subgroup, the same factor base, the same
targets and the same seed, and reports the total-operation ratio.  A
round whose lever is a cheaper restart, a cheaper probe or a cheaper
table leaves the variant's name alone, so that is the only pairing that
can see it.

The within-run pairing below is the other one: it prices a *new* row
against the rung it was built on, inside a single run.

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
across = json.load(open(arg("--across"))) if arg("--across") else None
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
            # A row whose column part repeats an earlier row's pins the
            # logarithm the way a rho collision does, not the way a
            # relation does.  Must be zero on every row of an honest
            # relation search; see the note's walk-merge diagnostic.
            ("repeated_column_rows", sum(r["linear_algebra"]["native"].get("repeated_column_rows", 0) for r in runs)),
            ("pinned_by_repeated_row", sum(r["linear_algebra"]["native"].get("pinned_by_repeated_row", 0) for r in runs)),
        ]))
    return out

def reproduction_check(new, old):
    """First-round rows of the new run against the frozen Round-1 file.

    A row whose counts are identical reproduces the earlier measurement,
    which is what makes the pairing below a comparison rather than two
    runs.  A row that moved is one the target guard bit: Round 1 drew
    targets without it, so a run could decompose one group element twice
    and be pinned by that collision instead of by its relations, ending
    early.  Those rows are reported with the size of the move, because
    the correction is to Round 1's numbers and not to this round's.
    """
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
            row = OrderedDict([
                ("regime", inst["regime"]),
                ("instance", inst["curve"]["name"]),
                ("variant", name),
                ("counts_identical", same),
            ])
            if not same:
                row["trials_round1"] = round(mean([r["trials"] for r in o]), 1)
                row["trials_guarded"] = round(mean([r["trials"] for r in runs]), 1)
                row["S_round1"] = round(mean([r["s"] for r in o]), 4)
                row["S_guarded"] = round(mean([r["s"] for r in runs]), 4)
                row["S_guarded_over_round1"] = round(mean([r["s"] for r in runs]) / mean([r["s"] for r in o]), 4)
                row["repeated_targets_skipped"] = round(mean([r["relations"]["native"].get("repeated_targets_skipped", 0) for r in runs]), 1)
            checks.append(row)
    return checks

def compare_across_rounds(new, old):
    """The same row, the same instance, the previous round's code.

    Pairs on `(regime, curve, variant, repeat index)`.  The two runs use
    the same seed, so repeat `i` of a row saw the same target in both,
    and the ratio of their total operation counts is a paired
    measurement rather than two independent draws.  A row present only
    in the candidate is reported as new — the within-run pairing above
    is what prices it.
    """
    out, new_rows = [], []
    old_by = {(i["regime"], i["curve"]["name"]): rows_of(i) for i in old["ledger"]["instances"]}
    for inst in new["ledger"]["instances"]:
        prev_rows = old_by.get((inst["regime"], inst["curve"]["name"]))
        if prev_rows is None:
            continue
        for name, runs in rows_of(inst).items():
            if name not in prev_rows:
                new_rows.append(OrderedDict([
                    ("regime", inst["regime"]),
                    ("instance", inst["curve"]["name"]),
                    ("variant", name),
                    ("S", round(mean([r["s"] for r in runs]), 4)),
                ]))
                continue
            base_runs = prev_rows[name]
            k = min(len(runs), len(base_runs))
            paired = [base_runs[i]["total_gae"] / runs[i]["total_gae"] for i in range(k)]
            # The phase the round's lever lives in, so a ratio that moves
            # can be read against the phase that moved it.
            phase = lambda rs, p: mean([r[p]["gae"] for r in rs])
            out.append(OrderedDict([
                ("regime", inst["regime"]),
                ("instance", inst["curve"]["name"]),
                ("log2_r", round(inst["log2_r"], 2)),
                ("variant", name),
                ("m", runs[0]["summands"]),
                ("signed_points", runs[0]["signed_points"]),
                ("columns", runs[0]["columns"]),
                ("base_changed", runs[0]["signed_points"] != base_runs[0]["signed_points"]),
                ("repeats_paired", k),
                ("speedup_per_repeat", [round(x, 4) for x in paired]),
                ("speedup_mean", round(mean(paired), 4)),
                ("speedup_min", round(min(paired), 4)),
                ("speedup_max", round(max(paired), 4)),
                ("S_previous_round", round(mean([r["s"] for r in base_runs]), 4)),
                ("S_this_round", round(mean([r["s"] for r in runs]), 4)),
                ("S_over_rho_previous_round", round(mean([r["s"] for r in base_runs]) / inst["rho_s_mean"], 4)),
                ("S_over_rho_this_round", round(mean([r["s"] for r in runs]) / inst["rho_s_mean"], 4)),
                ("relations_gae_previous_round", round(phase(base_runs, "relations"), 1)),
                ("relations_gae_this_round", round(phase(runs, "relations"), 1)),
                ("trials_previous_round", round(mean([r["trials"] for r in base_runs]), 1)),
                ("trials_this_round", round(mean([r["trials"] for r in runs]), 1)),
                ("walk_restarts_this_round", round(mean([r["relations"]["native"].get("walk_restarts", 0) for r in runs]), 1)),
                ("walk_jumps_this_round", round(mean([r["relations"]["native"].get("walk_jumps", 0) for r in runs]), 1)),
                ("walk_jumps_previous_round", round(mean([r["relations"]["native"].get("walk_jumps", 0) for r in base_runs]), 1)),
                ("ratio_to_family_optimum", round(mean([r.get("ratio_to_family_optimum", float("nan")) for r in runs]), 4)),
                ("verified_this_round", all(r["verified"] for r in runs)),
                ("verified_previous_round", all(r["verified"] for r in base_runs)),
                ("repeated_column_rows", sum(r["linear_algebra"]["native"].get("repeated_column_rows", 0) for r in runs)),
            ]))
    return out, new_rows

comparison = []
for inst in run["ledger"]["instances"]:
    comparison.extend(compare_instance(inst))

repro = reproduction_check(run, baseline) if baseline else []
across_rows, across_new = compare_across_rounds(run, across) if across else ([], [])
holdout_rows = []
if holdout:
    for inst in holdout["ledger"]["instances"]:
        holdout_rows.extend(compare_instance(inst))

# ── Markdown ──
if across_rows:
    print("**Against the previous round, same row, same instance, same seed.**\n")
    print("| regime | instance | log2 r | variant | m | |F| | speedup per repeat | speedup mean | min | max | S before → after | vs rho after | relations before → after | trials before → after | ok |")
    print("|:--|:--|--:|:--|--:|--:|:--|--:|--:|--:|:--|--:|:--|:--|:--|")
    for c in across_rows:
        print(
            f"| {c['regime']} | {c['instance']} | {c['log2_r']:.1f} | {c['variant']} | {c['m']} | {c['signed_points']:,} "
            f"| {', '.join(f(x) for x in c['speedup_per_repeat'])} | {f(c['speedup_mean'])} | {f(c['speedup_min'])} | {f(c['speedup_max'])} "
            f"| {f(c['S_previous_round'])} → {f(c['S_this_round'])} | {f(c['S_over_rho_this_round'])}× "
            f"| {f(c['relations_gae_previous_round'])} → {f(c['relations_gae_this_round'])} "
            f"| {c['trials_previous_round']:,.0f} → {c['trials_this_round']:,.0f} "
            f"| {'✓' if c['verified_this_round'] and c['verified_previous_round'] else '✗'} |"
        )
    if across_new:
        print(f"\nRows new in this round (priced by the within-run pairing below): "
              f"{', '.join(sorted({r['variant'] for r in across_new}))}.")
    print()
print("| regime | instance | log2 r | candidate | baseline | lever | class | m | |F| | K | speedup per repeat | speedup mean | vs first-round best | S before → after | vs rho after | vs floor after | y/c exact | ok |")
print("|:--|:--|--:|:--|:--|:--|:--|--:|--:|--:|:--|--:|--:|:--|--:|--:|--:|:--|")
for c in comparison:
    print(f"| {c['regime']} | {c['instance']} | {c['log2_r']:.1f} | {c['candidate']} | {c['baseline']} | {c['lever']} | {c['class']} | {c['m']} | {c['signed_points']} | {c['columns']} | {', '.join(f(x) for x in c['speedup_per_repeat'])} | {f(c['speedup_mean'])} | {f(c['speedup_vs_first_round_best'])} | {f(c['S_baseline'])} → {f(c['S_candidate'])} | {f(c['S_over_rho_candidate'])}× | {f(c['S_over_floor_candidate'])}× | {f(c['yield_over_ceiling_exact_candidate'])} | {'✓' if c['verified_candidate'] and c['verified_baseline'] and c['rho_verified'] else '✗'} |")
if repro:
    bad = [r for r in repro if not r["counts_identical"]]
    print(f"\nFirst-round rows against the frozen Round-1 file: {len(repro) - len(bad)} of {len(repro)} identical in every native count; {len(bad)} moved under the target guard.\n")
    if bad:
        print("| regime | instance | variant | trials Round 1 | trials guarded | S Round 1 | S guarded | ratio | repeats skipped |")
        print("|:--|:--|:--|--:|--:|--:|--:|--:|--:|")
        for r in bad:
            print(f"| {r['regime']} | {r['instance']} | {r['variant']} | {r['trials_round1']:,.0f} | {r['trials_guarded']:,.0f} | {f(r['S_round1'])} | {f(r['S_guarded'])} | {f(r['S_guarded_over_round1'])} | {r['repeated_targets_skipped']:.1f} |")
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
        ("previous_round", arg("--across")),
        ("holdout", arg("--holdout")),
        ("regression_suite_note", "The frozen WDSat solver regression suite (research/index_calculus_baseline_20260914/regression/) measures a SAT-solver stage on 60 fixed inputs and is not an adapter for other pipelines; the levers here (pair-table folds, target generation, exact ceiling, base sizing) touch no solver stage, so the matched full-DLP comparison in this file is the §8 evidence, under the parent accounting contract's exclusive per-phase charging."),
        ("all_candidates_verified", all(c["verified_candidate"] for c in comparison)),
        ("all_baselines_verified", all(c["verified_baseline"] for c in comparison)),
        ("no_row_pinned_by_a_repeated_column_part", all(c["repeated_column_rows"] == 0 for c in comparison + holdout_rows)),
        ("first_round_rows_identical", sum(1 for r in repro if r["counts_identical"]) if repro else None),
        ("first_round_rows_moved_by_the_target_guard", sum(1 for r in repro if not r["counts_identical"]) if repro else None),
        ("first_round_rows_moved_note", "Round 1 drew targets without the repeat guard, so a run could decompose one group element twice and be pinned by that collision rather than by its relations; the rows below are the ones that did, with the size of the correction to Round 1's figures"),
        ("reproduction_checks", repro),
        ("across_rounds_note", "Each row paired with the same-named row of the previous round's run on the same instance and seed; speedup = previous_round_total_operations / this_round_total_operations. This is the pairing that sees a lever which makes an existing row cheaper without renaming it."),
        ("across_rounds_all_verified", all(c["verified_this_round"] and c["verified_previous_round"] for c in across_rows) if across_rows else None),
        ("across_rounds", across_rows),
        ("across_rounds_new_variants", across_new),
        ("comparison", comparison),
        ("holdout_comparison", holdout_rows),
    ])
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print(f"\nwritten {out_path}")
