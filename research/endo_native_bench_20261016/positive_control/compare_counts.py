"""Compare the generic harness's CryptoPro-B counting record with the
specialised harness's committed one, arm by arm (stdlib only).

    python3 research/endo_native_bench_20261016/positive_control/compare_counts.py
"""
import json
import sys

A = "research/cryptopro_b_chain_sweep_20261006/counts.json"
B = "research/endo_native_bench_20261016/positive_control/cryptopro_b.counts.json"
FIELDS = ["mul", "sqr", "inv", "add", "m_eq", "m_eq_min", "m_eq_max", "correct"]

a, b = json.load(open(A)), json.load(open(B))
ka, kb = [x["key"] for x in a["arms"]], [x["key"] for x in b["arms"]]
print(f"arms: specialised {len(ka)}, generic {len(kb)}, same keys in the same order: {ka == kb}")
for k in ["chains_sha256", "pairs_m", "seed", "grid_chain", "element_chains", "reference",
          "continuity_chain", "continuity_config", "widths"]:
    print(f"inputs.{k} equal: {a['inputs'][k] == b['inputs'][k]}")
bb = {x["key"]: x for x in b["arms"]}
diffs = [(x["key"], f, x[f], bb[x["key"]][f]) for x in a["arms"] for f in FIELDS
         if x["key"] in bb and x[f] != bb[x["key"]][f]]
for d in diffs:
    print("DIFF", *d)
print(f"field differences over {len(ka)} arms x {len(FIELDS)} fields: {len(diffs)}")
print(f"all_correct: specialised {a['all_correct']}, generic {b['all_correct']}")
sys.exit(1 if diffs or ka != kb else 0)
