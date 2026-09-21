#!/usr/bin/env python3
"""Render a boundary-ledger JSON into the Markdown tables the note uses."""
import json, math, sys
from collections import OrderedDict

path = sys.argv[1]
d = json.load(open(path))
L = d["ledger"]

def f(v, digits=2):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "—"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    if abs(v) >= 100:
        return f"{v:.0f}"
    if abs(v) >= 10:
        return f"{v:.1f}"
    return f"{v:.{digits}f}"

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else float("nan")

print("## Instances\n")
print("| regime | instance | log2 r | #E/r | A | S_floor | rho S (mean) | rho walk S | rho steps/expected | rho ok |")
print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|")
for inst in L["instances"]:
    A = inst["automorphisms_generic"]
    ref = [x for x in inst["rho"] if (A <= 2 or x["automorphisms"] == A)]
    print(f"| {inst['regime']} | {inst['curve']['name']} | {inst['log2_r']:.1f} | {inst['cofactor']} | {A} | {inst['floor_s']:.3f} | {f(inst['rho_s_mean'])} | {f(mean([x['s_walk'] for x in ref]))} | {f(mean([x['steps_over_expected'] for x in ref]))} | {'✓' if inst['rho_verified_all'] else '✗'} |")

print("\n## Variants (means over repeats)\n")
print("| regime | instance | log2 r | variant | m | \\|F\\| | K | trials | yield/ceiling | S | S_wall | vs rho | vs floor | FB | rel | LA | verify | ok |")
print("|:--|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
for inst in L["instances"]:
    sqrt_r = math.sqrt(inst["r"])
    seen = OrderedDict()
    for v in inst["variants"]:
        seen.setdefault(v["name"], []).append(v)
    for name, runs in seen.items():
        m = lambda key: mean([r[key] for r in runs])
        ph = lambda p: mean([r[p]["gae"] for r in runs]) / sqrt_r
        ok = all(r["verified"] for r in runs)
        v0 = runs[0]
        print(f"| {inst['regime']} | {inst['curve']['name']} | {inst['log2_r']:.1f} | {name} | {v0['summands']} | {v0['signed_points']} | {v0['columns']} | {m('trials'):.0f} | {m('yield_over_ceiling'):.2f} | {f(m('s'))} | {f(m('s_wall'))} | {f(m('s')/inst['rho_s_mean'])} | {f(m('s')/inst['floor_s'])} | {f(ph('factor_base'))} | {f(ph('relations'))} | {f(ph('linear_algebra'))} | {f(ph('verify'), 3)} | {'✓' if ok else '✗'} |")

print("\n## Fits\n")
print("| regime | variant | phase | alpha (r) | R2 | alpha (#E) | R2 | points |")
print("|:--|:--|:--|--:|--:|--:|--:|--:|")
for fit in L["fits"]:
    print(f"| {fit['regime']} | {fit['variant']} | {fit['phase']} | {fit['alpha']:.3f} | {fit['r_squared']:.3f} | {fit['alpha_group_order']:.3f} | {fit['r_squared_group_order']:.3f} | {fit['points']} |")

op = d.get("oracle_pricing")
if op:
    print("\n## Oracle cells\n")
    print("| n | dim | m | \\|F\\| | K | unknowns | eqs | deg | eq/var | FFD | Macaulay D=2 rows×cols (rank) | D=3 | hit rate | disagreements |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|:--|--:|--:|")
    for c in op["cells"]:
        mac = {mm["degree"]: mm for mm in c["macaulay"]}
        cell = lambda D: f"{mac[D]['rows']}×{mac[D]['cols']} ({mac[D]['rank']})" if D in mac else "—"
        ffd = f"{c['ffd_min']}–{c['ffd_max']}" if c['ffd_min'] != c['ffd_max'] else f"{c['ffd_min']}"
        print(f"| {c['n']} | {c['dimension']} | {c['m']} | {c['signed_points']} | {c['signed_orbits']} | {c['unknowns']} | {c['equations']} | {c['system_degree']} | {c['eq_var_ratio']:.2f} | {ffd} | {cell(2)} | {cell(3)} | {c['hit_rate']:.2f} | {c['disagreements']} |")
    print("\n| n | m | oracle | unit | found/refuted/inconclusive | native median (found) | native median (refuted) | ms median (found) | ms median (refuted) | GAE/target mean | GAE/refutation mean | projected relation-phase S | vs floor |")
    print("|--:|--:|:--|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for c in op["cells"]:
        for o, p in zip(c["oracles"], c["projected"]):
            refut = f"{o['gae_per_refutation_mean']:.3e}" if o['refuted'] else "—"
            print(f"| {c['n']} | {c['m']} | {o['oracle']} | {o['native_unit']} | {o['found']}/{o['refuted']}/{o['inconclusive']} | {f(o['native_median_found'])} | {f(o['native_median_refuted'])} | {f(o['ms_median_found'],3)} | {f(o['ms_median_refuted'],3)} | {o['gae_per_target_mean']:.3e} | {refut} | {p['s_projected']:.3e} | {p['s_over_rho_floor']:.3e} |")
    ex = {}
    for c in op["cells"]:
        for o in c["oracles"]:
            if o["extra_totals"]:
                ex[(c['n'], c['m'], o['oracle'])] = o["extra_totals"]
    print("\nEngine detail totals:")
    for k, v in ex.items():
        print(" ", k, v)

print("\n## Calibration (ns per unit)\n")
print("| regime | instance | add | double | sqrt | AS solve | S4 pair | lookup | row op | word xor | legendre | inversion | frobenius |")
print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
for inst in L["instances"]:
    c = inst["calibration"]
    g = lambda k: f(c.get(k), 2) if c.get(k) is not None else "—"
    print(f"| {inst['regime']} | {inst['curve']['name']} | {g('ns_per_add')} | {g('ns_per_double')} | {g('ns_per_sqrt')} | {g('ns_per_as_solve')} | {g('ns_per_s4_pair')} | {g('ns_per_lookup')} | {g('ns_per_row_op')} | {g('ns_per_word_xor')} | {g('ns_per_legendre')} | {g('ns_per_inversion')} | {g('ns_per_frobenius')} |")

print("\nhost:", d.get("host"), "elapsed s:", d.get("elapsed_seconds"), "status:", d.get("status"), "all_verified:", d.get("all_verified"))
