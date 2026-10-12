#!/usr/bin/env python3
"""Registration R2 of RESEARCH_CACHEGRIND_N41_20260929.md: 15 randomized blocks of
{rho, IC} x {4 KiB, THP}, user/sys recorded separately, chase control before and after.
usage: run_r2.py         (env MEM_CALIB=<mem_calib binary>)     writes r2.json
"""
import json
import os
import random
import sys

import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("a2", os.path.join(HERE, "explore_stages_a2.py"))
a2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a2)
ex, ra = a2.ex, a2.ra
SEED, BLOCKS, MIB = 20260930, 15, 1024 * 1024


def one(arm, mode, tag):
    env_extra = ex.TUN if mode == "thp" else None
    ra.prepare()
    if arm == "rho":
        env = ra.rho_env()
        if env_extra:
            env.update(env_extra)
        rc, wall, u, s, huge = a2.run_polled2(ra.rho_cmd(), f"{ra.WORK}/rho_{tag}.jsonl", f"{ra.WORK}/rho_{tag}.stderr", env)
        return {"rc": rc, "wall_s": wall, "user_s": u, "sys_s": s, "max_anon_huge_kb": huge,
                "gates_G1": ra.rho_ok(f"{ra.WORK}/rho_{tag}.jsonl", rc)[0]}
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    rc, wall, u, s, huge = a2.run_polled2(ra.ic_cmd(f"{ra.WORK}/ic_{tag}.jsonl"), f"{ra.WORK}/ic_{tag}.summary.json", f"{ra.WORK}/ic_{tag}.stderr", env)
    ok, sm = ra.ic_ok(f"{ra.WORK}/ic_{tag}.summary.json", rc)
    return {"rc": rc, "wall_s": wall, "user_s": u, "sys_s": s, "max_anon_huge_kb": huge, "gates_G1": ok,
            "peak_rss_bytes": sm.get("peak_rss_bytes") if isinstance(sm, dict) else None}


def chase_pairs(n):
    return [{"k4": ex.chase(256 * MIB, 1, 4_000_000, 3), "thp": ex.chase(256 * MIB, 1, 4_000_000, 3, ex.TUN)} for _ in range(n)]


if __name__ == "__main__":
    if not ra.MEM:
        sys.exit("set MEM_CALIB")
    rng = random.Random(SEED)
    conds = [(a, m) for a in ("rho", "ic") for m in ("k4", "thp")]
    out = {"stage": "r2", "seed": SEED, "load_before": ra.uptime(), "chase_pre": chase_pairs(3), "blocks": []}
    for b in range(BLOCKS):
        order = conds[:]
        rng.shuffle(order)
        blk = {"order": [f"{a}/{m}" for a, m in order]}
        for a, m in order:
            blk[f"{a}/{m}"] = one(a, m, f"r2_b{b}_{a}_{m}")
        out["blocks"].append(blk)
        print(b, blk["order"], {k: round(v["user_s"], 2) for k, v in blk.items() if k != "order"}, flush=True)
    out["chase_post"] = chase_pairs(3)
    out["load_after"] = ra.uptime()
    json.dump(out, open(os.path.join(HERE, "r2.json"), "w"), indent=2)
