"""Render the lit2026 speedup journey from recorded entries.

Layout only: every plotted number is copied from an entry file, which copies it
from that step's run report. Nothing is computed here except axis placement.

usage: python -I render_journey.py <program_dir>
writes <program_dir>/journey/journey.svg and <program_dir>/journey/board.svg
"""
import json, pathlib, sys, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = pathlib.Path(sys.argv[1])
ws = json.loads((root / "workstreams.json").read_text())
entries = []
for p in sorted((root / "journey" / "entries").glob("*.json")):
    e = json.loads(p.read_text())
    e["_file"] = p.name
    entries.append(e)
superseded = {e.get("supersedes") for e in entries if e.get("supersedes")}
entries = [e for e in entries if e["_file"] not in superseded]

# ---- journey: per-workload speedup over the program baseline, by date ----
by_wl = collections.defaultdict(list)
for e in entries:
    r = (e.get("ratio_original_over_candidate") or {})
    if r.get("value") is not None:
        by_wl[(e["workload"]["id"], e["metric"], e["unit"])].append(e)

fig, ax = plt.subplots(figsize=(10, 5.4))
if by_wl:
    for (wid, metric, unit), es in sorted(by_wl.items()):
        es.sort(key=lambda e: (e["date"], e["_file"]))
        xs = list(range(len(es)))
        dates = [e["date"] for e in es]
        ys = [e["ratio_original_over_candidate"]["value"] for e in es]
        lo = [e["ratio_original_over_candidate"].get("interval", [None, None]) or [None, None] for e in es]
        line, = ax.plot(dates, ys, marker="o", label=f"{wid} ({metric}, {unit})")
        for d, y, (a, b) in zip(dates, ys, lo):
            if a is not None and b is not None:
                ax.vlines(d, a, b, color=line.get_color(), alpha=0.5)
        for d, y, e in zip(dates, ys, es):
            ax.annotate(f'{e["workstream"]} {e["slug"]}', (d, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.axhline(1.0, color="#555", lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel("program baseline / candidate (same session; >1 is faster)")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=7)
else:
    ax.text(0.5, 0.5, "No measured steps recorded yet.\nEach workstream PR adds journey/entries/*.json\nwith its interleaved original/previous/candidate arms.",
            ha="center", va="center", transform=ax.transAxes, fontsize=11, color="#555")
    ax.set_xticks([]); ax.set_yticks([])
n_all = len(entries)
ax.set_title(f"lit2026 speedup journey: measured speedup over program baseline {ws['program_baseline_commit']}\n"
             f"({n_all} recorded step(s); every ratio measured in one interleaved session with the baseline binary)", fontsize=10)
fig.text(0.01, 0.01, "Source: research/lit2026_speedup_program/journey/entries/*.json", fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(root / "journey" / "journey.svg")

# ---- board: workstream status and recorded steps by kind ----
kinds = ["speedup", "no-change", "regression", "diagnostic", "negative"]
kcol = {"speedup": "#3aa76d", "no-change": "#9aa0a6", "regression": "#d0453b", "diagnostic": "#2a6fdb", "negative": "#7a4cc2"}
cnt = collections.Counter((e["workstream"], e["kind"]) for e in entries)
ids = [w["id"] for w in ws["workstreams"]]
fig, ax = plt.subplots(figsize=(10, 5.2))
left = [0] * len(ids)
for k in kinds:
    vals = [cnt[(i, k)] for i in ids]
    ax.barh(ids, vals, left=left, color=kcol[k], label=k, edgecolor="white")
    left = [a + b for a, b in zip(left, vals)]
for i, w in enumerate(ws["workstreams"]):
    ax.text(max(left) + 0.3 if max(left) else 0.1, i, f'[{w["status"]}] {w["title"][:78]}', va="center", fontsize=6.5)
ax.invert_yaxis()
ax.set_xlim(0, (max(left) if max(left) else 0) + 9)
ax.set_xlabel("recorded journey steps (count)")
ax.set_title("lit2026 workstreams: status and recorded steps by kind", fontsize=10)
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.xaxis.get_major_locator().set_params(integer=True)
fig.text(0.01, 0.01, "Source: research/lit2026_speedup_program/workstreams.json and journey/entries/*.json", fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(root / "journey" / "board.svg")
print(f"rendered {n_all} entries")
