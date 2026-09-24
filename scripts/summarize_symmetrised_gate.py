#!/usr/bin/env python3
"""X4' of research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTE_TARGETS.md, applied as
registered to the frozen output of examples/koblitz_symmetrised_gate.rs.

Every cost is in group-addition equivalents (GAE) per relation found.  An
oracle's cost is its word XORs (elimination and specialisation, a lower
bound) priced at the raw word-XOR rate over one curve addition, both
measured at the rung.  The reference is the cheaper enumeration rule on the
same base: `full` (every pair, every decomposition harvested) or `first_hit`
(stop at the first, pay the full enumeration where there is none).

The gate is closed if the symmetrised oracle costs at least the reference at
n = 23 and the slope of log2(Q_sym / Q_ref) against n is not negative with a
two-standard-error interval excluding zero.  The frozen calibration
(docs/ic/calibration.json, K_0 entries only) re-prices the word XOR where it
has an entry; if that changes the verdict, the gate is undetermined by the
conversion.  The 350x question is the slope of log2(Q_sym / Q_x), per
relation, on the d = 3 diagonal, with the wall-clock ratio's slope beside it;
a rung where either arm is budget-limited on more than half its targets is
printed as a bound and left out of that fit.

Several gate files may be given (the d = 3 run and the separately run,
non-deciding d = 4 rungs); a JSON object with `instances` is read as the
calibration instead of docs/ic/calibration.json.  The verdicts use d = 3
only, as registered.

Usage:
    python3 scripts/summarize_symmetrised_gate.py GATE.json [GATE.json ...] [CALIBRATION.json]
"""

from __future__ import annotations

import json
import math
import sys

E1 = {13, 19, 23}


def fit(xs, ys):
    """Least-squares slope and its standard error."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    if n < 3:
        return b, float("nan")
    res = [y - (my + b * (x - mx)) for x, y in zip(xs, ys)]
    return b, math.sqrt(sum(e * e for e in res) / (n - 2) / sxx)


def reference(e, gae_per_step):
    full = e["full_steps"] * e["_targets"] * gae_per_step / e["decompositions_total"] if e["decompositions_total"] else math.inf
    first = e["first_hit_steps_total"] * gae_per_step / e["decomposable"] if e["decomposable"] else math.inf
    return full, first, min(full, first)


def arm_q(a, targets, gae_per_xor):
    total = a["mean_word_xors"] * targets * gae_per_xor
    return total / a["found"] if a["found"] else math.inf


def gate_verdict(ns, ratios):
    xs = [float(n) for n in ns]
    ys = [math.log2(r) for r in ratios]
    b, se = fit(xs, ys)
    top = ratios[ns.index(max(ns))]
    falling = b < 0 and not math.isnan(se) and b + 2 * se < 0
    return ("closed" if top >= 1 and not falling else "open"), b, se, top


def main() -> None:
    rows, calib = [], None
    for path in sys.argv[1:]:
        doc = json.load(open(path))
        if isinstance(doc, dict) and "instances" in doc:
            calib = doc["instances"]
        else:
            rows.extend(doc)
    if calib is None:
        calib = json.load(open("docs/ic/calibration.json"))["instances"]
    table = []
    for b in rows:
        k = b["targets"]
        enum = {e["base"]: dict(e, _targets=k) for e in b["enumeration"]}
        cal = calib.get(f"koblitz/K_{b['a']} / GF(2^{b['n']})")
        by_cap = {}
        for a in b["arms"]:
            by_cap.setdefault(a["cap"] if a["arm"] == "x-chained" else a["cap"] - 1, {})[a["arm"]] = a
        for d, arms in sorted(by_cap.items()):
            x, s = arms["x-chained"], arms["symmetrised"]
            fx, hx, rx = reference(enum["x"], enum["x"]["gae_per_step"])
            fu, hu, ru = reference(enum["u"], enum["u"]["gae_per_step"])
            qx, qs = arm_q(x, k, b["gae_per_word_xor"]), arm_q(s, k, b["gae_per_word_xor"])
            row = dict(n=b["n"], l=b["l"], d=d, x=x, s=s, fx=fx, hx=hx, rx=rx, fu=fu, hu=hu, ru=ru, qx=qx, qs=qs,
                       gae_xor=b["gae_per_word_xor"], gae_xor_rref=b["ns_per_word_xor_in_rref"] / b["ns_per_add"],
                       fx_pt=x["mean_word_xors"] * b["gae_per_word_xor"] / (enum["x"]["full_steps"] * enum["x"]["gae_per_step"]),
                       fu_pt=s["mean_word_xors"] * b["gae_per_word_xor"] / (enum["u"]["full_steps"] * enum["u"]["gae_per_step"]),
                       wall_x=x["mean_ms"] * k / x["found"] if x["found"] else math.inf,
                       wall_s=s["mean_ms"] * k / s["found"] if s["found"] else math.inf,
                       censored=x["budget"] > k / 2 or s["budget"] > k / 2,
                       gates=x["gate_failures"] + s["gate_failures"], cal=cal, k=k,
                       enum_u=enum["u"], enum_x=enum["x"])
            if cal:
                # The frozen ratio re-prices the word XOR; a step is one addition and one lookup.
                step = 1.0 + cal["ns_per_lookup"]
                _, _, ru_c = reference(enum["u"], step)
                row["ratio_cal"] = arm_q(s, k, cal["ns_per_word_xor"]) / ru_c
            table.append(row)

    print("| n | shape | l | d | arm | vars | top-degree only | found / refuted / budget | gate failures "
          "| GAE per relation | reference (full, first-hit) | ratio to reference | per target vs full | splits | ms per relation "
          "| built degree | oversize targets |")
    print("|---:|---|---:|---:|---|---:|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|")
    for r in table:
        for name, a, q, ref, full, first, pt, wall in [
            ("x-chained", r["x"], r["qx"], r["rx"], r["fx"], r["hx"], r["fx_pt"], r["wall_x"]),
            ("symmetrised", r["s"], r["qs"], r["ru"], r["fu"], r["hu"], r["fu_pt"], r["wall_s"]),
        ]:
            shape = "E1" if r["n"] in E1 else "fit only"
            print(f"| {r['n']} | {shape} | {r['l']} | {r['d']} | {name} (cap {a['cap']}) | {a['n_vars']} "
                  f"| {'yes' if a['top_degree_only'] else 'no'} | {a['found']} / {a['refuted']} / {a['budget']} "
                  f"| {a['gate_failures']} | {q:,.0f} | {ref:,.0f} ({full:,.0f}, {first:,.0f}) | {q / ref:.2f} "
                  f"| {pt:.2f} | {a['mean_splits']:.0f} | {wall:.1f} | {a['built_degree']} | {a['oversize_targets']} |")
    print()
    for r in table:
        cal = f", frozen {r['cal']['ns_per_word_xor']:.4f}" if r["cal"] else ", no frozen K_0 entry"
        print(f"n = {r['n']}: {r['gae_xor']:.4f} GAE per word XOR at the raw rate (in the solver's elimination "
              f"{r['gae_xor_rref']:.4f}{cal}); |F_x| = {r['enum_x']['points']}, |F_u| = {r['enum_u']['points']}; "
              f"decomposable {r['enum_x']['decomposable']}/{r['k']} on x, {r['enum_u']['decomposable']}/{r['k']} on u")
    if any(r["gates"] for r in table):
        print("\nGATE FAILURES: the rungs above with a non-zero count are invalid.")

    d3 = [r for r in table if r["d"] == 3 and not r["gates"]]
    print()
    ns = [r["n"] for r in d3]
    ratios = [r["qs"] / r["ru"] for r in d3]
    verdict, b, se, top = gate_verdict(ns, ratios)
    print(f"gate (d = 3, {len(ns)} rungs): Q_sym / Q_ref = " + ", ".join(f"{n}: {q:.2f}" for n, q in zip(ns, ratios)))
    print(f"  slope of log2(Q_sym/Q_ref) in n: {b:+.4f} ± {se:.4f} (1 s.e.); at n = {max(ns)}: {top:.2f}  ->  {verdict.upper()}")
    if any("ratio_cal" in r for r in d3):
        alt = [r.get("ratio_cal", q) for r, q in zip(d3, ratios)]
        v2, b2, se2, top2 = gate_verdict(ns, alt)
        print(f"  with the frozen calibration where it has an entry: " + ", ".join(f"{n}: {q:.2f}" for n, q in zip(ns, alt))
              + f"; slope {b2:+.4f} ± {se2:.4f}  ->  {v2.upper()}")
        if v2 != verdict:
            print("  UNDETERMINED BY THE CONVERSION")

    fitted = [r for r in d3 if not r["censored"]]
    bounds = [r for r in d3 if r["censored"]]
    xs = [float(r["n"]) for r in fitted]
    b, se = fit(xs, [math.log2(r["qs"] / r["qx"]) for r in fitted])
    bw, sew = fit(xs, [math.log2(r["wall_s"] / r["wall_x"]) for r in fitted])
    print()
    print("350x question (d = 3): Q_sym / Q_x per relation = "
          + ", ".join(f"{r['n']}: {r['qs'] / r['qx']:.2f}" for r in fitted)
          + ("; bounds (a budget-limited arm): " + ", ".join(f"{r['n']}: {r['qs'] / r['qx']:.2f}" for r in bounds) if bounds else ""))
    print(f"  slope in n: {b:+.4f} ± {se:.4f}; wall-clock ratio slope {bw:+.4f} ± {sew:.4f}")
    if len(fitted) < 3 or math.isnan(se):
        cls = "too few unbounded rungs to fit"
    elif b < 0 and b + 2 * se < 0:
        cls = "ADVANCE CANDIDATE" if bw < 0 else "UNDETERMINED (the two partial units disagree)"
    else:
        cls = "ENGINEERING (as predicted)"
    print(f"  -> {cls}")


if __name__ == "__main__":
    main()
