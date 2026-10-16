"""Tabulate the committed jacobian-transfers.json against its previous version.

Extraction and counting only: every class below was computed by the native
`jacobian_transfers build`; this script reads two committed files and counts.

usage: python3 tabulate.py OLD_JSON NEW_JSON OUT_DIR
Writes route_changes.csv (cover_descent rows repriced at the GHS floor, by k)
and verdict_changes.csv (curve verdicts, old -> new).
"""
import collections
import csv
import hashlib
import json
import sys

old_path, new_path, out = sys.argv[1:4]
old, new = json.load(open(old_path)), json.load(open(new_path))
assert [c["slug"] for c in old["curves"]] == [c["slug"] for c in new["curves"]]

routes, verdicts = collections.Counter(), collections.Counter()
for a, b in zip(old["curves"], new["curves"]):
    verdicts[(a["verdict"], b["verdict"])] += 1
    for r, s in zip(a["routes"], b["routes"]):
        inv = s["invariants"]
        if s["kind"] != "cover_descent" or "class_at_genus_k" not in inv:
            continue
        # the old file and the new file's own record of the old class agree
        assert r["class"] == inv["class_at_genus_k"], a["slug"]
        k = a["field"]["degree"]
        routes[(k, r["destination"]["genus"], s["destination"]["genus"], r["class"], s["class"])] += 1

with open(f"{out}/route_changes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["k", "genus_before", "genus_after", "class_before", "class_after", "routes"])
    for key, n in sorted(routes.items()):
        w.writerow([*key, n])
with open(f"{out}/verdict_changes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["verdict_before", "verdict_after", "curves"])
    for key, n in sorted(verdicts.items()):
        w.writerow([*key, n])
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()  # noqa: E731
print("old", sha(old_path)[:16], "new", sha(new_path)[:16])
