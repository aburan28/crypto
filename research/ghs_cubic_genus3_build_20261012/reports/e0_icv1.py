"""Compute the ICV1 identity of each source E0 model from the run records.
Independent order check: scripts/curve_id.py counts #E(GF(p^3)) itself (Euler's criterion)
and must agree with order_E0 from the Rust tool. Derivation only.
Usage (repo root): python3 research/ghs_cubic_genus3_build_20261012/reports/e0_icv1.py"""
import json, sys
sys.path.insert(0, "scripts")
import curve_id as cid
D = "research/ghs_cubic_genus3_build_20261012/results/"
out = []
for run in ["run-20261013T0320Z", "run-20261013T0347Z"]:
    for line in open(D + run + ".jsonl"):
        r = json.loads(line)
        if r["rec"] != "case":
            continue
        w = r["E0_weierstrass_ext"]
        order = cid._count_extension(w["p"], w["modulus"], w["a"], w["b"])
        assert order == w["order"] == r["order_E0"], (r["u_idx"], order)
        ident = cid.extension_id(w["p"], 3, w["modulus"], w["a"], w["b"], order)
        out.append({"run": run, "q": r["q"], "u_idx": r["u_idx"], "kind": r["kind"],
                    "field_modulus_low_coefficients": w["modulus"], "a": w["a"], "b": w["b"],
                    "order": order, "icv1": ident["icv1"], "slug": ident["slug"],
                    "model_json": ident["model_json"]})
json.dump({"what": "Weierstrass model Y^2=X^3-27IX-27J of the Jacobian of the quartic y0^2=u*u'' (isomorphic to E0 over K, since the quartic has K-points at infinity); order independently recounted by scripts/curve_id.py",
           "registered_in_docs_curves_registry": False, "curves": out},
          open(D + "e0_icv1.json", "w"), indent=1)
for o in out:
    print(o["run"], o["q"], o["u_idx"], o["slug"])
