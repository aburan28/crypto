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
}
CHIP = {"advance": "advance, count", "relabel": "relabelling"}
GROUP = {"prime": "Prime field, generic curve", "char2": "Binary field, random curve", "koblitz": "Koblitz curve, signed-Frobenius reference"}

by_regime = OrderedDict()
for inst in L["instances"]:
    by_regime.setdefault(inst["regime"], []).append(inst)

print("<!-- table rows: largest instance per regime -->")
for regime, insts in by_regime.items():
    big = max(insts, key=lambda i: i["r"])
    print(f'          <tr class="group"><td colspan="6">{GROUP[regime]} &middot; {html.escape(big["curve"]["name"])} &middot; r = 2<sup>{big["log2_r"]:.1f}</sup> &middot; rho S = {fmt(big["rho_s_mean"])} &middot; floor S = {big["floor_s"]:.3f}</td></tr>')
    seen = OrderedDict()
    for v in big["variants"]:
        seen.setdefault(v["name"], []).append(v)
    for name, runs in seen.items():
        s = mean([r["s"] for r in runs])
        cls, why = CLASS.get((regime, name), CLASS.get(name, ("baseline", "baseline")))
        chip = f'<span class="chip {cls}">{CHIP[cls]}</span>' if cls in CHIP else f'<span class="chip">{cls}</span>'
        ok = "yes" if all(r["verified"] for r in runs) else "no"
        print(f'          <tr><td>{PRETTY.get(name, name)}</td><td class="n">{fmt(s)}</td><td class="n">{fmt(s / big["rho_s_mean"])}&times;</td><td class="n">{fmt(s / big["floor_s"])}&times;</td><td>{ok}</td><td title="{html.escape(why)}">{chip}</td></tr>')

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
