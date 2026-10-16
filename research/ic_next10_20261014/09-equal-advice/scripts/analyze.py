"""Item 9 analysis: equal-advice IC vs BL rho, plain rho reference, census.

Run: python3 -I scripts/analyze.py   (from the item directory)
Writes results.json, fig_online_vs_advice.svg/.png, fig_census.svg/.png.
"""
import glob
import json
import math
import os
import random
import statistics as st

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
B = 10000
random.seed(20261014)


def jl(path):
    return [json.loads(l) for l in open(path) if l.strip()]


def ic_bytes(K, n):
    cand = K * K * n
    slots = 1 << (4 * cand - 1).bit_length()
    return 16 * slots + 24 * cand + 32 * 2 * K * n + 8 * K * n


def boot_ratio(num, den):
    """Paired target bootstrap: ratio of means and median of per-target ratios."""
    k = len(num)
    rm, md = [], []
    for _ in range(B):
        idx = [random.randrange(k) for _ in range(k)]
        a = [num[i] for i in idx]
        b = [den[i] for i in idx]
        rm.append(sum(a) / sum(b))
        md.append(st.median([x / y for x, y in zip(a, b)]))
    rm.sort()
    md.sort()
    lo, hi = int(0.025 * B), int(0.975 * B) - 1
    return {
        "ratio_of_means": sum(num) / sum(den),
        "ratio_of_means_ci95": [rm[lo], rm[hi]],
        "median_ratio": st.median([x / y for x, y in zip(num, den)]),
        "median_ratio_ci95": [md[lo], md[hi]],
    }


def load_ic(name):
    rows = [d for d in jl(os.path.join(DATA, name)) if d.get("kind") == "compact_orbit_dlp_target"]
    rows.sort(key=lambda d: d["fixture_index"])
    return rows


def load_bl(name):
    rows = jl(os.path.join(DATA, name))
    pre = [d for d in rows if d["kind"] == "bl_precompute"][0]
    levels = {}
    for d in rows:
        if d["kind"] == "bl_online_target":
            levels.setdefault(d["prefix_level"], []).append(d)
    for v in levels.values():
        v.sort(key=lambda d: d["target_index"])
    return pre, levels


def load_rho(n):
    rows = []
    for path in sorted(glob.glob(os.path.join(DATA, f"rho_n{n}_t*.json"))):
        recs = [d for d in jl(path) if d.get("kind") == "rho_ks_batch_fixture"]
        if recs:
            rows.append(recs[0])
    return rows


def scalars(n):
    return [int(l) for l in open(os.path.join(DATA, f"scalars_n{n}.txt")) if l.strip()]


def summarize(xs):
    return {"mean": st.mean(xs), "median": st.median(xs), "min": min(xs), "max": max(xs), "count": len(xs)}


results = {"rungs": {}}
plan = [
    (41, [100], "bl_n41.jsonl"),
    (53, [200], "bl_n53.jsonl"),
    (61, [256, 360], "bl_n61_64th.jsonl"),
]
for n, ks, blname in plan:
    rung = {}
    truth = scalars(n)
    ic_runs = {}
    for K in ks:
        path = os.path.join(DATA, f"ic_n{n}_K{K}.jsonl")
        if not os.path.exists(path):
            continue
        rows = load_ic(f"ic_n{n}_K{K}.jsonl")
        if len(rows) == 0:
            continue
        ok = sum(1 for d, t in zip(rows, truth) if d["group_verified"] and d["recovered_scalar"] == t)
        ic_runs[K] = rows
        rung[f"ic_K{K}"] = {
            "advice_bytes": ic_bytes(K, n),
            "targets": len(rows),
            "verified_vs_truth": ok,
            "online_ms": summarize([d["online_ms"] for d in rows]),
            "probes": summarize([d["probes"] for d in rows]),
            "root_calls": summarize([d["root_calls"] for d in rows]),
            "ns_per_root_call": 1e6 * sum(d["online_ms"] for d in rows) / sum(d["root_calls"] for d in rows),
        }
    blpath = os.path.join(DATA, blname)
    if os.path.exists(blpath):
        pre, levels = load_bl(blname)
        rung["bl_precompute"] = pre
        for lv, rows in sorted(levels.items()):
            solved = [d for d in rows if d["solved"]]
            ok = sum(1 for d in rows if d["solved"] and d["recovered_scalar"] == truth[d["target_index"]])
            entry = {
                "table_entries": rows[0]["table_entries"],
                "table_bytes": rows[0]["table_bytes"],
                "targets": len(rows),
                "solved": len(solved),
                "verified_vs_truth": ok,
                "online_ms": summarize([d["online_ms"] for d in rows]),
                "online_steps": summarize([d["online_steps"] for d in rows]),
                "walks": summarize([d["walks"] for d in rows]),
                "ns_per_step": 1e6 * sum(d["online_ms"] for d in rows) / max(1, sum(d["online_steps"] for d in rows)),
                "predicted_steps_sqrt_N_over_T": math.sqrt(pre["classes_N"] / rows[0]["table_entries"]),
            }
            for K, icrows in ic_runs.items():
                if len(icrows) == len(rows):
                    entry[f"ic_K{K}_over_bl"] = boot_ratio(
                        [d["online_ms"] for d in icrows], [d["online_ms"] for d in rows]
                    )
                    entry[f"ic_K{K}_bytes_over_bl_bytes"] = ic_bytes(K, n) / rows[0]["table_bytes"]
            rung[f"bl_level{lv}"] = entry
        # positive control: log-log slope of mean steps vs entries
        pts = [(math.log(v[0]["table_entries"]), math.log(st.mean(d["online_steps"] for d in v))) for v in levels.values()]
        if len(pts) >= 2:
            mx = st.mean(p[0] for p in pts)
            my = st.mean(p[1] for p in pts)
            slope = sum((x - mx) * (y - my) for x, y in pts) / sum((x - mx) ** 2 for x, _ in pts)
            rung["bl_positive_control_slope"] = slope
    mm = os.path.join(DATA, f"bl_n{n}_mismatch.jsonl")
    if os.path.exists(mm):
        pre_m, lv_m = load_bl(f"bl_n{n}_mismatch.jsonl")
        rows = lv_m.get(0, [])
        rung["bl_negative_control_mismatched_map"] = {
            "targets": len(rows),
            "solved": sum(1 for d in rows if d["solved"]),
            "online_steps": summarize([d["online_steps"] for d in rows]) if rows else None,
            "budget_steps": rows[0]["budget_steps"] if rows else None,
        }
    rho = load_rho(n)
    if rho:
        k = len(rho)
        rung["plain_strong_rho"] = {
            "targets": k,
            "verified": sum(1 for d in rho if d.get("verified")),
            "online_ms": summarize([d["online_ms"] for d in rho]),
            "walk_steps": summarize([d["walk_steps"] for d in rho]),
        }
        for K, icrows in ic_runs.items():
            rung["plain_strong_rho"][f"rho_over_ic_K{K}"] = boot_ratio(
                [d["online_ms"] for d in rho], [d["online_ms"] for d in icrows[:k]]
            )
        if "bl_level0" in rung:
            blrows = load_bl(blname)[1][0]
            rung["plain_strong_rho"]["rho_over_bl_level0"] = boot_ratio(
                [d["online_ms"] for d in rho], [d["online_ms"] for d in blrows[:k]]
            )
    results["rungs"][str(n)] = rung

# census
cen = []
for path in sorted(glob.glob(os.path.join(DATA, "census_n*_K*.jsonl"))):
    s = [d for d in jl(path) if d["kind"] == "free_relation_census"]
    if s:
        cen.append(s[0])
results["census"] = cen
json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=1)

# ---------- figure 1: online ms vs advice bytes
fig, ax = plt.subplots(figsize=(8, 5))
colors = {41: "#4c78a8", 53: "#f58518", 61: "#54a24b"}
for n, ks, _ in plan:
    rung = results["rungs"].get(str(n), {})
    c = colors[n]
    bl = sorted([(v["table_bytes"], v["online_ms"]["mean"]) for k, v in rung.items() if k.startswith("bl_level")])
    if bl:
        ax.plot([b / 2**20 for b, _ in bl], [m for _, m in bl], "o-", color=c, label=f"n={n} BL preprocessing rho")
    for K in ks:
        e = rung.get(f"ic_K{K}")
        if e:
            ax.plot(e["advice_bytes"] / 2**20, e["online_ms"]["mean"], "D" if K == 360 else "s", color=c, markersize=10,
                    markeredgecolor="black", label=f"n={n} IC K={K}")
    r = rung.get("plain_strong_rho")
    if r:
        ax.axhline(r["online_ms"]["mean"], color=c, ls=":", lw=1)
        ax.text(0.6, r["online_ms"]["mean"] * 1.08, f"n={n} plain rho (no advice)", color=c, fontsize=8)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("advice (MiB)")
ax.set_ylabel("mean online wall per target (ms), 64 targets")
ax.set_title("Online cost against advice: compact-orbit IC vs Bernstein-Lange rho")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=7, loc="lower left")
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(os.path.join(HERE, f"fig_online_vs_advice.{ext}"), dpi=150)

# ---------- figure 2: census
if cen:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for d in cen:
        c = colors.get(d["n"], "gray")
        ax.plot(d["expected_pairs_uniform"], max(d["verified_nontrivial"], 0.5), "o", color=c)
        ax.annotate(f"n={d['n']} K={d['K']}", (d["expected_pairs_uniform"], max(d["verified_nontrivial"], 0.5)),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.plot(d["leads_prediction_n4K4_over_r"], max(d["verified_nontrivial"], 0.5), "x", color=c, alpha=0.5)
        ax.plot(d["expected_pairs_uniform"], max(d["rank_mod_r"], 0.5), "^", color=c, alpha=0.8)
    lim = [0.5, max(max(d["leads_prediction_n4K4_over_r"] for d in cen), 200) * 1.5]
    ax.plot(lim, lim, "k--", lw=1, label="measured = predicted")
    ax.axhline(100, color="red", lw=1, ls=":", label="kill line: 100 relations at n=61")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("predicted pairs (o: uniform M^2 n / r; x: leads note n^4K^4/r)")
    ax.set_ylabel("verified non-trivial collisions (o, x); rank mod r (^)")
    ax.set_title("E-3 free-relation census")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(HERE, f"fig_census.{ext}"), dpi=150)
print(json.dumps(results, indent=1)[:6000])
