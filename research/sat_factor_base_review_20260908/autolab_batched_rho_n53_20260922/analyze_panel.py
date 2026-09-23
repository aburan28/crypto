#!/usr/bin/env python3
import json, re, statistics as st, pathlib, sys
LAB = pathlib.Path(__file__).resolve().parent
D = LAB / "runs" / "panel_v2"
def wall(p): return float(re.search(r"([\d.]+) real", p.read_text()).group(1))
def rss(p): return int(re.search(r"(\d+)\s+maximum resident", p.read_text()).group(1))
def lines(p, kind): return [json.loads(l) for l in p.open() if f'"kind":"{kind}"' in l]
out = {"panel": {}, "complete": True}
for L in (1024, 2048, 4096, 8192):
    rows = []
    for b in (0, 1, 2):
        icp, ksp = D / f"ic{L}_b{b}", D / f"ks{L}_b{b}"
        if not (icp.with_suffix(".time").exists() and ksp.with_suffix(".time").exists()):
            out["complete"] = False; continue
        ic = lines(icp.with_suffix(".jsonl"), "relation_rank_summary")
        ks = lines(ksp.with_suffix(".jsonl"), "rho_ks_batch_fixture")
        ks_sum = lines(ksp.with_suffix(".jsonl"), "rho_ks_batch_summary")
        ic_sum = lines(icp.with_suffix(".jsonl"), "retained_support_batch_summary")
        match = len(ic) == len(ks) == L and all(
            i["published_fixture_scalar"] == k["published_fixture_scalar"] == i["recovered_fixture_scalar"] == k["recovered_fixture_scalar"]
            and i["published_q"] == k["published_q"] for i, k in zip(ic, ks))
        icw, ksw = wall(icp.with_suffix(".time")), wall(ksp.with_suffix(".time"))
        rows.append(dict(block=b, order="ic_first" if b % 2 == 0 else "rho_first", targets_match=match,
            ic_all_solved=bool(ic_sum and ic_sum[0]["all_fixtures_solved"]),
            ic_wall_s=icw, rho_wall_s=ksw, ic_over_rho=icw / ksw,
            ic_peak_rss=rss(icp.with_suffix(".time")), rho_peak_rss=rss(ksp.with_suffix(".time")),
            ic_charged_ms=ic_sum[0]["full_algorithm_charged_total_ms"] if ic_sum else None,
            rho_group_additions=ks_sum[0]["charges"]["group_additions"] if ks_sum else None,
            load_ic=(icp.with_suffix(".load")).read_text().strip(), load_rho=(ksp.with_suffix(".load")).read_text().strip()))
    if not rows: continue
    r = [x["ic_over_rho"] for x in rows]
    verdict = "ic_overtakes" if len(rows) == 3 and all(x < 1 for x in r) else ("rho_faster" if len(rows) == 3 and all(x > 1 for x in r) else "unresolved")
    out["panel"][L] = dict(rows=rows, median_ic_over_rho=st.median(r), min=min(r), max=max(r), verdict=verdict,
        all_match=all(x["targets_match"] for x in rows),
        rho_additions_deterministic=len({x["rho_group_additions"] for x in rows}) == 1)
    print(f"L={L:5d} blocks={len(rows)} IC/rho median={st.median(r):.3f} [{min(r):.3f},{max(r):.3f}] "
          f"IC wall med={st.median(x['ic_wall_s'] for x in rows):.1f}s rho={st.median(x['rho_wall_s'] for x in rows):.1f}s "
          f"match={out['panel'][L]['all_match']} verdict={verdict}")
(LAB / "runs" / "panel_v2_summary.json").write_text(json.dumps(out, indent=1))
