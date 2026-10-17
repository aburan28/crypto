"""Post-hoc (not preregistered): compare recorded dispersion indices with the exact binomial null 1 - q,
q = per-x probability of landing in a given box; and count corner-density intervals excluding 1. Reads results.json only."""
import json, sys
d = sys.argv[1] if len(sys.argv) > 1 else "."
r = json.load(open(f"{d}/results.json")); p = r["pair46"]["p"]
out = {"note": "post hoc; exact null for box counts is Binomial(B1, q), var/mean = 1 - q", "curves": []}
for c in r["pair46"]["curves"]:
    rows = []
    for dd in c["dispersion"]:
        q = 0.5 * dd["H"] / ((p + 1) // 2)  # P(liftable) * P(cy in band)
        null = 1 - q
        rows.append({"B1": dd["B1"], "index": dd["index"], "ci95": dd["ci95"], "exact_null": null,
                     "contains_exact_null": dd["ci95"][0] <= null <= dd["ci95"][1],
                     "index_over_null": dd["index"] / null})
    ex = [e["k"] for e in c["corner"] if e["k"] > 0 and not (e["ci95"][0] <= 1 <= e["ci95"][1])]
    out["curves"].append({"name": c["name"], "dispersion": rows, "corner_k_excluding_1": ex})
json.dump(out, open(f"{d}/posthoc.json", "w"), indent=1)
for c in out["curves"]:
    print(c["name"], [(x["B1"], round(x["index_over_null"], 4), x["contains_exact_null"]) for x in c["dispersion"]], "corner k excl 1:", c["corner_k_excluding_1"])
