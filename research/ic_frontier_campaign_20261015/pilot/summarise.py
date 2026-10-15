"""Lay out the Phase 0 pilot records as a Markdown table (recorded values only).

    python3 summarise.py sessions/*/records.jsonl > TABLE.md
"""
import json
import sys

print("| session | curve | log2 r | arm | status | S | S complete? | unpriced counters (count or 'not emitted') |")
print("|:--|:--|--:|:--|:--|--:|:--|:--|")
for path in sys.argv[1:]:
    for line in open(path):
        r = json.loads(line)
        if r.get("warmup"):
            continue
        c = r["cost"]
        w = r["workload"]["curve"]
        import math
        lr = math.log2(int(w["r"]))
        un = c.get("unpriced") or []
        shown = []
        for k in un:
            v = r["counters"].get(k)
            shown.append(f"{k.replace('_uncharged', '')}={'n/e' if v is None else v}")
        s = c.get("s")
        print(f"| {r['session_id']} | {w['slug']} | {lr:.2f} | {r['arm']} | {r['outcome']['status']} | "
              f"{'' if s is None else f'{s:.3f}'} | {'no (lower bound)' if c['lower_bound'] else 'yes'} | "
              f"{', '.join(shown) if shown else '—'} |")
