#!/usr/bin/env python3
"""Render a round's tables for the boundary-ledger note.

usage: boundary_round_note_tables.py <round.json> [--ladder]

The first table is one block per regime: the generic floor, the family
optimum, the reference, every first-round row of the largest instance,
and every later rung with the `S` of the rung it was built on as its
"was" mark.  The second is the fitted total exponent of every new
variant next to the row it was built on and the reference.  `--ladder`
adds the full ladder of new rows.
"""
import json, math, sys
from collections import OrderedDict

run = json.load(open(sys.argv[1]))
L = run["ledger"]
ROUND2 = ("_negfold", "_frobfold", "_walk", "_balanced")
SUFFIXES = ("_balanced", "_walk", "_frobfold", "_negfold")

PRETTY = {
    "semaev_s3_roots_m2": "Semaev `S₃` roots",
    "direct_subtraction_m2": "direct subtraction",
    "mitm_m2": "meet in the middle",
    "mitm_m3": "meet in the middle",
    "semaev_s4_pairs_and_solve_m3": "`S₄` pairs-and-solve",
    "mitm_m2_negfold": "+ negation-folded table",
    "mitm_m2_negfold_walk": "+ walk targets",
    "mitm_m3_negfold": "+ negation-folded table",
    "mitm_m3_negfold_walk": "+ walk targets",
    "mitm_m3_signed_orbit_columns": "meet in the middle, signed-orbit columns",
    "mitm_m2_signed_orbit_columns": "meet in the middle, signed-orbit columns",
    "mitm_m3_abscissa_columns_control": "the same base, one column per abscissa",
    "mitm_m2_abscissa_columns_control": "the same base, one column per abscissa",
    "semaev_s4_pairs_and_solve_m3_signed_orbit_columns": "`S₄` pairs-and-solve on the invariant subspace",
    "mitm_m3_signed_orbit_columns_negfold": "+ negation-folded table",
    "mitm_m2_signed_orbit_columns_negfold": "+ negation-folded table",
    "mitm_m3_signed_orbit_columns_frobfold": "+ Frobenius-folded table",
    "mitm_m2_signed_orbit_columns_frobfold": "+ Frobenius-folded table",
    "mitm_m3_signed_orbit_columns_frobfold_walk": "+ walk targets",
    "mitm_m2_signed_orbit_columns_frobfold_walk": "two summands, folded table, walk targets",
    "mitm_m3_signed_orbit_columns_frobfold_walk_balanced": "balanced base, three summands",
    "mitm_m2_signed_orbit_columns_frobfold_walk_balanced": "balanced base, two summands",
    "mitm_m2_negfold_walk_balanced": "+ balanced base",
}
CLASS = {
    "semaev_s3_roots_m2": "baseline",
    "direct_subtraction_m2": "accounting",
    "mitm_m2": "engineering",
    "mitm_m3": "engineering",
    "semaev_s4_pairs_and_solve_m3": "relabelling",
    "mitm_m3_signed_orbit_columns": "advance, count",
    "mitm_m2_signed_orbit_columns": "advance, count",
    "mitm_m3_abscissa_columns_control": "baseline (control)",
    "mitm_m2_abscissa_columns_control": "baseline (control)",
    "semaev_s4_pairs_and_solve_m3_signed_orbit_columns": "relabelling",
}
LABEL = {"prime": "**prime**", "char2": "**binary**", "koblitz": "**Koblitz**"}

def f(v, d=2):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "—"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    if abs(v) >= 100:
        return f"{v:.0f}"
    if abs(v) >= 10:
        return f"{v:.1f}"
    return f"{v:.{d}f}"

def mean(xs):
    return sum(xs) / len(xs)

def rows_of(inst):
    seen = OrderedDict()
    for v in inst["variants"]:
        seen.setdefault(v["name"], []).append(v)
    return seen

def previous_rung(name, names):
    for s in SUFFIXES:
        if name.endswith(s) and name[: -len(s)] in names:
            return name[: -len(s)]
    return None

by_regime = OrderedDict()
for inst in L["instances"]:
    by_regime.setdefault(inst["regime"], []).append(inst)

print("| regime, instance | variant | m | \\|F\\| | K | trials | yield/ceiling (exact) | S | was | vs rho | vs floor | vs family | ok | class |")
print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|")
for regime, insts in by_regime.items():
    big = max(insts, key=lambda i: i["r"])
    A = big["automorphisms_generic"]
    ref = [x for x in big["rho"] if (A <= 2 or x["automorphisms"] == A)]
    walk_s = mean([x["s_walk"] for x in ref])
    head = f"{LABEL[regime]}, `{big['curve']['name']}`, `r = 2^{big['log2_r']:.1f}`, `#E = {big['cofactor']}r`, `A = {A}`"
    fam = big["variants"][0].get("family_optimum_s")
    print(f"| {head} | generic floor `√(π/2A)` | | | | | | {big['floor_s']:.3f} | | {f(big['floor_s'] / big['rho_s_mean'])}× | 1× | | — | boundary |")
    if fam:
        print(f"| | family optimum `0.75·#E^(2/3)/(a√r)` at `\\|F\\| = #E^(1/3)` | | {big['variants'][0]['family_optimum_base']:,.0f} | | | | {f(fam)} | | {f(fam / big['rho_s_mean'])}× | {f(fam / big['floor_s'])}× | 1× | — | model |")
    name_ref = "signed-Frobenius rho, counted" if regime == "koblitz" else "Pollard rho, r-adding, counted"
    print(f"| | {name_ref} (walk alone {f(walk_s)}) | | | | | | {f(big['rho_s_mean'])} | | 1× | {f(big['rho_s_mean'] / big['floor_s'])}× | | ✓ | reference |")
    rows = rows_of(big)
    first = [n for n in rows if not any(s in n for s in ROUND2)]
    best_before = min(first, key=lambda n: mean([r["s"] for r in rows[n]]))
    best_after = min((n for n in rows if any(s in n for s in ROUND2)),
                     key=lambda n: mean([r["s"] for r in rows[n]]), default=None)
    for name, runs in rows.items():
        s = mean([r["s"] for r in runs])
        is_r2 = any(x in name for x in ROUND2)
        was = ""
        if is_r2:
            p = previous_rung(name, rows) or best_before
            was = f(mean([r["s"] for r in rows[p]]))
        cls = "engineering" if is_r2 else CLASS.get(name, "baseline")
        mark = (lambda t: f"**{t}**") if name in (best_before, best_after) else (lambda t: t)
        print(
            f"| | {PRETTY.get(name, name)} | {runs[0]['summands']} | {runs[0]['signed_points']:,} | {runs[0]['columns']} "
            f"| {mean([r['trials'] for r in runs]):,.0f} | {mean([r['yield_over_ceiling'] for r in runs]):.2f} "
            f"({mean([r['yield_over_ceiling_exact'] for r in runs]):.2f}) | {mark(f(s))} | {was} "
            f"| {mark(f(s / big['rho_s_mean']) + '×')} | {f(s / big['floor_s'])}× "
            f"| {f(mean([r.get('ratio_to_family_optimum', float('nan')) for r in runs]))}× "
            f"| {'✓' if all(r['verified'] for r in runs) else '✗'} | {cls} |"
        )

print()
print("| regime | variant | α (r) | R² | sizes | α of the rung it was built on |")
print("|:--|:--|--:|--:|--:|--:|")
fits = {(x["regime"], x["variant"], x["phase"]): x for x in L["fits"]}
for regime, insts in by_regime.items():
    big = max(insts, key=lambda i: i["r"])
    rows = rows_of(big)
    first = [n for n in rows if not any(s in n for s in ROUND2)]
    best_before = min(first, key=lambda n: mean([r["s"] for r in rows[n]]))
    rho = fits.get((regime, "rho_reference", "total"))
    if rho:
        print(f"| {regime} | rho reference | {rho['alpha']:.3f} | {rho['r_squared']:.3f} | {rho['points']} | — |")
    for name in rows:
        fit = fits.get((regime, name, "total"))
        if not fit:
            continue
        p = previous_rung(name, rows) or (best_before if any(s in name for s in ROUND2) else None)
        pf = fits.get((regime, p, "total")) if p else None
        print(
            f"| {regime} | {name} | {fit['alpha']:.3f} | {fit['r_squared']:.3f} | {fit['points']} "
            f"| {f(pf['alpha'], 3) if pf else '—'} |"
        )

if "--ladder" in sys.argv:
    print()
    print("| regime | instance | log₂ r | variant | m | \\|F\\| | K | table | targets | trials | y/c | exact | S | was | vs rho | vs floor | FB | rel | LA | ok |")
    print("|:--|:--|--:|:--|--:|--:|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for inst in L["instances"]:
        sq = math.sqrt(inst["r"])
        rows = rows_of(inst)
        first = [n for n in rows if not any(s in n for s in ROUND2)]
        best_before = min(first, key=lambda n: mean([r["s"] for r in rows[n]]))
        for name, runs in rows.items():
            if not any(x in name for x in ROUND2):
                continue
            s = mean([r["s"] for r in runs])
            p = previous_rung(name, rows) or best_before
            ph = lambda k: mean([r[k]["gae"] for r in runs]) / sq
            print(
                f"| {inst['regime']} | {inst['curve']['name']} | {inst['log2_r']:.1f} | {name} | {runs[0]['summands']} "
                f"| {runs[0]['signed_points']} | {runs[0]['columns']} | {runs[0]['table'].replace('_folded', '')} | {runs[0]['targets']} "
                f"| {mean([r['trials'] for r in runs]):,.0f} | {mean([r['yield_over_ceiling'] for r in runs]):.2f} "
                f"| {mean([r['yield_over_ceiling_exact'] for r in runs]):.2f} | {f(s)} | {f(mean([r['s'] for r in rows[p]]))} "
                f"| {f(s / inst['rho_s_mean'])} | {f(s / inst['floor_s'])} | {f(ph('factor_base'))} | {f(ph('relations'))} "
                f"| {f(ph('linear_algebra'))} | {'✓' if all(r['verified'] for r in runs) else '✗'} |"
            )
