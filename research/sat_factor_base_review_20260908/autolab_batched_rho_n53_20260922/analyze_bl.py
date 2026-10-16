#!/usr/bin/env python3
import json, re, statistics as st, pathlib
LAB = pathlib.Path(__file__).resolve().parent; D = LAB / "runs" / "bl_warm"
def q(xs, p): xs = sorted(xs); return xs[int(p * (len(xs) - 1))]
def rss(p): return int(re.search(r"(\d+)\s+maximum resident", p.read_text()).group(1))
def wall(p): return float(re.search(r"([\d.]+) real", p.read_text()).group(1))
ic = [json.loads(l) for l in open(D / "ic.jsonl") if '"relation_rank_summary"' in l]
ic_pre = (ic[0]["curve_setup_ms"] + ic[0]["setup_ms"] + ic[0]["collection_ms"] + ic[0]["linear_solve_ms"]) / 1e3
ic_warm = [f["collection_ms"] + f["linear_solve_ms"] for f in ic[1:]]
out = {"ic": dict(precompute_s=ic_pre, warm_median_ms=st.median(ic_warm), warm_mean_ms=st.mean(ic_warm),
       p95=q(ic_warm, .95), p99=q(ic_warm, .99), max=max(ic_warm), peak_rss=rss(D / "ic.time"), wall_s=wall(D / "ic.time"))}
print(f"IC  precompute {ic_pre:6.1f}s  warm median {out['ic']['warm_median_ms']:6.2f} mean {out['ic']['warm_mean_ms']:6.2f} p95 {out['ic']['p95']:6.2f} p99 {out['ic']['p99']:6.2f} max {out['ic']['max']:7.2f} ms  RSS {out['ic']['peak_rss']/1e9:.2f} GB")
for P in (242000, 969000, 3875000):
    f = D / f"bl_{P}.jsonl"
    if not f.exists() or not (D / f"bl_{P}.time").exists(): continue
    rows = [json.loads(l) for l in open(f)]
    fx = [r for r in rows if r.get("kind") == "rho_ks_batch_fixture"]; s = rows[-1]
    if s.get("kind") != "rho_ks_batch_summary": continue
    ics = {tuple(i["published_q"]): i["published_fixture_scalar"] for i in ic}
    match = len(fx) == 1024 and all(ics.get(tuple(x["published_q"])) == x["published_fixture_scalar"] == x["recovered_fixture_scalar"] for x in fx)
    w = [x["total_ms"] for x in fx]
    pre = (s["setup_ms"] + s["precompute_ms"]) / 1e3
    out[P] = dict(precompute_s=pre, precompute_over_ic=pre / ic_pre, precompute_steps=s["precompute_steps"], table=s["precompute_table_entries"],
                  warm_median_ms=st.median(w), warm_mean_ms=st.mean(w), p95=q(w, .95), p99=q(w, .99), max=max(w),
                  all_via_precomputed=all(x["solved_via_precomputed"] for x in fx), targets_match_ic=match, peak_rss=rss(D / f"bl_{P}.time"))
    o = out[P]
    print(f"BL P={P:>8} precompute {pre:6.1f}s ({o['precompute_over_ic']:.2f}x IC)  warm median {o['warm_median_ms']:6.2f} mean {o['warm_mean_ms']:6.2f} p95 {o['p95']:6.2f} p99 {o['p99']:6.2f} max {o['max']:7.2f} ms  table {o['table']}  RSS {o['peak_rss']/1e6:.0f} MB  match={match}")
(LAB / "runs" / "bl_warm_summary.json").write_text(json.dumps(out, indent=1, default=str))
