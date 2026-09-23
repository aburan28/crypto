import json, re, glob, statistics as st, os
os.chdir(os.path.dirname(__file__) + "/runs/paired")
def wall(p): return float(re.search(r"([\d.]+) real", open(p).read()).group(1))
def rss(p): return int(re.search(r"(\d+)\s+maximum resident", open(p).read()).group(1))
def lines(p, kind): return [json.loads(l) for l in open(p) if f'"kind":"{kind}"' in l]
out = {}
for L in (32, 1024):
    rows = []
    for b in (0, 1, 2):
        ic = lines(f"ic{L}_b{b}.jsonl", "relation_rank_summary")
        ks = lines(f"ks{L}_b{b}.jsonl", "rho_ks_batch_fixture")
        assert len(ic) == len(ks) == L
        match = all(i["published_fixture_scalar"] == k["published_fixture_scalar"] == i["recovered_fixture_scalar"] == k["recovered_fixture_scalar"]
                    and i["published_q"] == k["published_q"] for i, k in zip(ic, ks))
        icw, ksw = wall(f"ic{L}_b{b}.time"), wall(f"ks{L}_b{b}.time")
        rows.append(dict(block=b, order="ic_first" if b % 2 == 0 else "rho_first", targets_match=match,
                         ic_wall_s=icw, ks_wall_s=ksw, ic_over_ks=icw / ksw,
                         ic_rss=rss(f"ic{L}_b{b}.time"), ks_rss=rss(f"ks{L}_b{b}.time"),
                         ks_steps=sum(k["walk_steps"] for k in ks)))
    out[L] = rows
    print(f"L={L}")
    for r in rows: print("  ", r)
    print("   median ic/ks wall ratio:", round(st.median(r["ic_over_ks"] for r in rows), 3),
          " all targets match:", all(r["targets_match"] for r in rows),
          " ks steps identical across blocks:", len({r["ks_steps"] for r in rows}) == 1)
json.dump(out, open("../paired_summary.json", "w"), indent=1)
