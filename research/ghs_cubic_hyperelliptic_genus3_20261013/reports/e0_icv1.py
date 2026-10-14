"""ICV1 identity of each E0 (short Weierstrass model over K in the census basis) from the run records.
Naming only, as in research/ghs_cubic_genus3_build_20261012/reports/e0_icv1.py: scripts/curve_id.py
recounts #E(GF(p^3)) independently and must agree with order_E0 from the Rust tool.
Usage (repo root): python3 research/ghs_cubic_hyperelliptic_genus3_20261013/reports/e0_icv1.py"""
import json, sys
sys.path.insert(0, "scripts")
import curve_id as cid
D = "research/ghs_cubic_hyperelliptic_genus3_20261013/results/"
RUN = "run-20261014T0029Z"
mods, out = {}, []
for line in open(D + RUN + ".jsonl"):
    r = json.loads(line)
    if r["rec"] == "field":
        mods[r["q"]] = r["census_K_modulus_low"]
    if r["rec"] != "instance":
        continue
    q, w = r["q"], r["short_weierstrass_census_basis"]
    order = cid._count_extension(q, mods[q], w["A"], w["B"])
    assert order == r["order_E0"], (q, r["Q"], order)
    ident = cid.extension_id(q, 3, mods[q], w["A"], w["B"], order)
    out.append({"run": RUN, "q": q, "arm": r["arm"], "Q": r["Q"], "p_census": r["p_census"],
                "field_modulus_low_coefficients": mods[q], "A": w["A"], "B": w["B"], "order": order,
                "icv1": ident["icv1"], "slug": ident["slug"]})
json.dump({"what": "short Weierstrass model y^2 = x^3 + A x + B over GF(q)[t]/(t^3 + c2 t^2 + c1 t + c0) of E0 (K-isomorphic to y0^2 = F0(x) via the root r0); order recounted by scripts/curve_id.py",
           "registered_in_docs_curves_registry": False, "curves": out}, open(D + "e0_icv1.json", "w"), indent=1)
for o in out:
    print(o["q"], o["arm"], o["Q"], o["order"], o["slug"])
