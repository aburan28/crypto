#!/usr/bin/env python3
"""Summarize result JSONL files: verified / attempts and total seconds."""
import collections, json, sys

cells = collections.OrderedDict()
for path in sys.argv[1:]:
    for line in open(path):
        r = json.loads(line)
        if "meta" in r:
            continue
        key = (r["n"], r["budget_s"], "planted" if r["planted"] else "random",
               r["formulation"])
        c = cells.setdefault(key, {"ok": 0, "k": 0, "t": 0.0, "unsat": 0, "rej": 0})
        c["k"] += 1
        c["ok"] += r["status"] == "verified"
        c["unsat"] += r["status"] == "unsat"
        c["rej"] += r.get("rejected_tuples", 0)
        c["t"] += r["total_s"]
print("| n | budget s | targets | formulation | verified | total s | unsat | rejected tuples |")
print("|---|---:|---|---|---:|---:|---:|---:|")
for (n, b, kind, form), c in cells.items():
    print(f"| {n} | {b:g} | {kind} | {form} | {c['ok']}/{c['k']} | {c['t']:.2f} | {c['unsat']} | {c['rej']} |")
