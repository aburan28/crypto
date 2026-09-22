#!/usr/bin/env python3
"""Emit the HTML fragments for the scoreboard's boundary-ledger panel."""
import json, math, sys, html
from collections import OrderedDict

run = json.load(open(sys.argv[1]))
L = run["ledger"]

def fmt(v):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "&mdash;"
    if v >= 1000:
        return f"{v:,.0f}"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 10:
        return f"{v:.1f}"
    return f"{v:.2f}"

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else float("nan")

PRETTY = {
    "semaev_s3_roots_m2": "Semaev S&#8323; roots, m = 2",
    "direct_subtraction_m2": "Direct subtraction, m = 2",
    "mitm_m2": "Meet in the middle, m = 2",
    "mitm_m3": "Meet in the middle, m = 3",
    "semaev_s4_pairs_and_solve_m3": "S&#8324; pairs-and-solve, m = 3",
    "mitm_m3_signed_orbit_columns": "Meet in the middle, signed-orbit columns",
    "mitm_m2_signed_orbit_columns": "Meet in the middle (m = 2), signed-orbit columns",
    "mitm_m3_abscissa_columns_control": "Same base, one column per abscissa (control)",
    "mitm_m2_abscissa_columns_control": "Same base (m = 2), one column per abscissa (control)",
    "semaev_s4_pairs_and_solve_m3_signed_orbit_columns": "S&#8324; pairs-and-solve on the invariant subspace",
    # Round 2: the cumulative engineering ledger, one rung per suffix.
    "mitm_m2_negfold": "Meet in the middle, m = 2, negation-folded table",
    "mitm_m2_negfold_walk": "&hellip; and walk targets",
    "mitm_m3_negfold": "Meet in the middle, m = 3, negation-folded table",
    "mitm_m3_negfold_walk": "&hellip; and walk targets",
    "mitm_m3_signed_orbit_columns_negfold": "Signed-orbit columns, negation-folded table",
    "mitm_m3_signed_orbit_columns_frobfold": "&hellip; Frobenius-folded table (|F|&sup2;/4n entries)",
    "mitm_m3_signed_orbit_columns_frobfold_walk": "&hellip; and walk targets",
    "mitm_m2_signed_orbit_columns_negfold": "Signed-orbit columns (m = 2), negation-folded table",
    "mitm_m2_signed_orbit_columns_frobfold": "&hellip; Frobenius-folded table (|F|&sup2;/4n entries)",
    "mitm_m2_signed_orbit_columns_frobfold_walk": "Two summands, Frobenius-folded table, walk targets",
    "mitm_m3_signed_orbit_columns_frobfold_walk_balanced": "Balanced base (|F| &asymp; 1.2&middot;(4#E)<sup>1/3</sup>), m = 3, folded table, walk",
    "mitm_m2_signed_orbit_columns_frobfold_walk_balanced": "Balanced base, m = 2, folded table, walk",
    # Round 3: the balanced base on the other two regimes, sized by the
    # family shape law rather than by the ceil(bits/3) rule.
    "mitm_m2_negfold_walk_balanced": "&hellip; and a base at the family optimum (|F| = #E<sup>1/3</sup>)",
}
# Class of each variant by the AGENTS.md §3 test, argued in research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md §3.
# Keyed by (regime, variant); a bare variant name is the fallback for every regime.
CLASS = {
    "semaev_s3_roots_m2": ("baseline", "baseline: one square root per base abscissa per target"),
    "direct_subtraction_m2": ("accounting", "accounting: the same enumeration, priced by additions instead of square roots"),
    "mitm_m2": ("engineering", "engineering: memory for time, ratio to the counting boundary unmoved"),
    "mitm_m3": ("engineering", "engineering: memory for time, one probe per subtraction"),
    ("char2", "mitm_m3"): ("baseline", "baseline: the table pipeline on the random binary curve"),
    "semaev_s4_pairs_and_solve_m3": ("relabel", "relabelling: the pair table's memory is paid back as a per-target loop over |F|^2/2 pairs; S rises"),
    "mitm_m3_signed_orbit_columns": ("advance", "advance, count: K falls by the orbit length 2n against the control; S moves little because the factor base dominates"),
    "mitm_m2_signed_orbit_columns": ("advance", "advance, count: K falls by the orbit length 2n against the control; S moves little because the factor base dominates"),
    "mitm_m3_abscissa_columns_control": ("baseline", "baseline: the random-binary pipeline on the Koblitz curve"),
    "mitm_m2_abscissa_columns_control": ("baseline", "baseline: the random-binary pipeline on the Koblitz curve"),
    "semaev_s4_pairs_and_solve_m3_signed_orbit_columns": ("relabel", "relabelling: pairs instead of probes on the same base"),
    # Round 2 (RESEARCH_IC_BOUNDARY_LEDGER.md §10): every rung is engineering by the §3 test —
    # S falls, trials, yield/ceiling and the exponents do not move.
    "mitm_m2_negfold": ("engineering", "engineering: half the table's additions build entries the abscissa key already held; same trials, same yield"),
    "mitm_m2_negfold_walk": ("engineering", "engineering: one addition per target instead of two scalar multiplications; same yield against the ceiling"),
    "mitm_m3_negfold": ("engineering", "engineering: the table's additions halve; same trials, same yield"),
    "mitm_m3_negfold_walk": ("engineering", "engineering: the target's two scalar multiplications become one addition; the |F| subtractions per target remain"),
    "mitm_m3_signed_orbit_columns_negfold": ("engineering", "engineering: the table's additions halve; same trials, same yield"),
    "mitm_m3_signed_orbit_columns_frobfold": ("engineering", "engineering: the table shrinks by the orbit length 2n, the automorphism group the floor already credits to a generic walk; one canonicalisation per probe, priced"),
    "mitm_m3_signed_orbit_columns_frobfold_walk": ("engineering", "engineering: one addition per target; at m = 3 the |F| subtractions per target dominate, so little moves"),
    "mitm_m2_signed_orbit_columns_negfold": ("engineering", "engineering: the table's additions halve; same trials, same yield"),
    "mitm_m2_signed_orbit_columns_frobfold": ("engineering", "engineering: the table shrinks by the orbit length 2n; one canonicalisation per probe, priced"),
    "mitm_m2_signed_orbit_columns_frobfold_walk": ("engineering", "engineering: two summands on the folded table with walk targets; the trials are what the exact ceiling allows, at one addition each"),
    "mitm_m3_signed_orbit_columns_frobfold_walk_balanced": ("engineering", "engineering: a different base, sized so the folded table balances the trials; its own row, its own floor"),
    "mitm_m2_signed_orbit_columns_frobfold_walk_balanced": ("engineering", "engineering: the same balanced base with two summands; the best Koblitz row, still above the reference"),
    # Round 3 (§11): still engineering — the base size is a parameter of the
    # same counting bound, now set by the family shape law instead of a rule.
    "mitm_m2_negfold_walk_balanced": ("engineering", "engineering: the base size that minimises the family's table-plus-relations cost, F = #E^(1/3) signed points; its own row, its own table"),
}
CHIP = {"advance": "advance, count", "relabel": "relabelling"}
ROUND2 = ("_negfold", "_frobfold", "_walk", "_balanced")

def previous_rung(name, names):
    """The rung this Round-2 row was built on: the name with its last suffix
    stripped, else the first-round best of the instance."""
    for s in ("_walk", "_frobfold", "_negfold", "_balanced"):
        if name.endswith(s) and name[: -len(s)] in names:
            return name[: -len(s)]
    return None
GROUP = {"prime": "Prime field, generic curve", "char2": "Binary field, random curve", "koblitz": "Koblitz curve, signed-Frobenius reference"}

by_regime = OrderedDict()
for inst in L["instances"]:
    by_regime.setdefault(inst["regime"], []).append(inst)

print("<!-- table rows: largest instance per regime -->")
for regime, insts in by_regime.items():
    big = max(insts, key=lambda i: i["r"])
    print(f'          <tr class="group"><td colspan="7">{GROUP[regime]} &middot; {html.escape(big["curve"]["name"])} &middot; r = 2<sup>{big["log2_r"]:.1f}</sup> &middot; rho S = {fmt(big["rho_s_mean"])} &middot; floor S = {big["floor_s"]:.3f}</td></tr>')
    seen = OrderedDict()
    for v in big["variants"]:
        seen.setdefault(v["name"], []).append(v)
    first_round = [n for n in seen if not any(s in n for s in ROUND2)]
    best_before = min(first_round, key=lambda n: mean([r["s"] for r in seen[n]])) if first_round else None
    for name, runs in seen.items():
        s = mean([r["s"] for r in runs])
        cls, why = CLASS.get((regime, name), CLASS.get(name, ("baseline", "baseline")))
        chip = f'<span class="chip {cls}">{CHIP[cls]}</span>' if cls in CHIP else f'<span class="chip">{cls}</span>'
        ok = "yes" if all(r["verified"] for r in runs) else "no"
        was = ""
        if any(sfx in name for sfx in ROUND2):
            prev = previous_rung(name, seen) or best_before
            if prev:
                was = f' <small>was {fmt(mean([r["s"] for r in seen[prev]]))}</small>'
        # The family shape law is a model of this family, not a bound on the
        # problem (§11.2), so it gets its own column and its own legend line
        # rather than sitting beside the floor as if it were one.
        fam = mean([r.get("ratio_to_family_optimum") for r in runs])
        print(f'          <tr><td>{PRETTY.get(name, name)}</td><td class="n">{fmt(s)}{was}</td><td class="n">{fmt(s / big["rho_s_mean"])}&times;</td><td class="n">{fmt(s / big["floor_s"])}&times;</td><td class="n">{fmt(fam)}&times;</td><td>{ok}</td><td title="{html.escape(why)}">{chip}</td></tr>')

print("\n<!-- exponent rows: total per variant, r-fit -->")
def pos(alpha):
    return max(0.0, min(100.0, (alpha - 0.25) / 0.65 * 100.0))
for f in L["fits"]:
    if f["phase"] != "total":
        continue
    name = f["variant"]
    label = f'{f["regime"]}: {PRETTY.get(name, name) if name != "rho_reference" else "rho reference"}'
    print(f'        <div class="erow"><div class="name">{label}</div><div class="track"><i class="dot" style="left:{pos(f["alpha"]):.1f}%" title="r^{f["alpha"]:.2f} (R² {f["r_squared"]:.2f}, {f["points"]} sizes)"></i></div><div class="num">{f["alpha"]:.2f}</div></div>')

print("\n<!-- phase exponents (relations / linear algebra), for the note -->")
for f in L["fits"]:
    if f["phase"] in ("relations", "linear_algebra", "factor_base"):
        print(f'{f["regime"]} {f["variant"]} {f["phase"]}: alpha_r={f["alpha"]:.3f} R2={f["r_squared"]:.3f} alpha_E={f["alpha_group_order"]:.3f} n={f["points"]}')
