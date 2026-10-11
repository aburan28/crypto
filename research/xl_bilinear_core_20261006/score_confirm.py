"""Score results/confirm.json against PROTOCOL.md X1, X2 (thresholds as frozen)."""
import json

res = json.load(open("results/confirm.json"))
cells = res["cells"]
off = [c for c in cells if c["bmax"] < c["l"]]
cap = [c for c in cells if c["bmax"] == c["l"]]
within = [c for c in off if abs(c["bmax"] - c["tight"]) <= 1]
cap_ok = [c for c in cap if c["tight"] >= c["l"] - 2]
edge = [c for c in cells if c["bmax"] in (c["tested"][0], c["tested"][-1]) and c["bmax"] < c["l"]
        and c["bmax"] != 0]
incons = sum(r["inconsistent"] for r in res["raw"])
score = dict(
    X1=dict(off_cap_cells=len(off), within1=len(within),
            frac=round(len(within) / len(off), 3) if off else None,
            misses=[(c["n"], c["l"], (c["D1"], c["D2"]), c["bmax"], c["tight"]) for c in off if c not in within],
            cap_cells=len(cap), cap_ok=len(cap_ok),
            cap_misses=[(c["n"], c["l"], (c["D1"], c["D2"]), c["bmax"], c["tight"]) for c in cap if c not in cap_ok],
            bmax_at_tested_edge=[(c["n"], c["l"], (c["D1"], c["D2"]), c["bmax"]) for c in edge]),
    X2=dict(inconsistent_planted=incons, pass_=incons == 0),
)
score["X1"]["pass_"] = bool(off) and len(within) / len(off) >= 0.9 and len(cap_ok) == len(cap)
json.dump(score, open("results/score_confirm.json", "w"), indent=1)
print(json.dumps(score, indent=1))
