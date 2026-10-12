#!/usr/bin/env python3
"""Amendment A2 of RESEARCH_CACHEGRIND_N41_20260929.md: the E2 huge-page test again, with user and
system CPU time recorded separately (A1's E2 recorded only their sum).

usage: explore_stages_a2.py           (env MEM_CALIB=<mem_calib binary>)
7 interleaved repetitions of {rho, IC} x {4 KiB pages, GLIBC_TUNABLES=glibc.malloc.hugetlb=1}, core 0.
Writes explore_thp_a2.json.
"""
import json
import os
import subprocess
import sys
import time

import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ex", os.path.join(HERE, "explore_stages.py"))
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)
ra = ex.ra


def run_polled2(cmd, out, err, env, pin=0):
    full = ["taskset", "-c", str(pin)] + cmd
    t0 = time.perf_counter()
    with open(out, "w") as fo, open(err, "w") as fe:
        p = subprocess.Popen(full, stdout=fo, stderr=fe, env=env)
        best = 0
        while True:
            done, status, ru = os.wait4(p.pid, os.WNOHANG)
            if done:
                break
            v = ex.anon_huge_kb(p.pid)
            if v:
                best = max(best, v)
            time.sleep(0.25)
    return os.waitstatus_to_exitcode(status), time.perf_counter() - t0, ru.ru_utime, ru.ru_stime, best


def pair(tag, env_extra=None):
    ra.prepare()
    rec = {}
    renv = ra.rho_env()
    ienv = dict(os.environ)
    if env_extra:
        renv.update(env_extra)
        ienv.update(env_extra)
    rc, wall, u, s, huge = run_polled2(ra.rho_cmd(), f"{ra.WORK}/rho_{tag}.jsonl", f"{ra.WORK}/rho_{tag}.stderr", renv)
    rec["rho"] = {"rc": rc, "wall_s": wall, "user_s": u, "sys_s": s, "gates_G1": ra.rho_ok(f"{ra.WORK}/rho_{tag}.jsonl", rc)[0], "max_anon_huge_kb": huge}
    rc, wall, u, s, huge = run_polled2(ra.ic_cmd(f"{ra.WORK}/ic_{tag}.jsonl"), f"{ra.WORK}/ic_{tag}.summary.json", f"{ra.WORK}/ic_{tag}.stderr", ienv)
    ok, sm = ra.ic_ok(f"{ra.WORK}/ic_{tag}.summary.json", rc)
    rec["ic"] = {"rc": rc, "wall_s": wall, "user_s": u, "sys_s": s, "gates_G1": ok, "max_anon_huge_kb": huge,
                 "peak_rss_bytes": sm.get("peak_rss_bytes") if isinstance(sm, dict) else None}
    return rec


if __name__ == "__main__":
    if not ra.MEM:
        sys.exit("set MEM_CALIB")
    out = {"stage": "thp_a2", "load_before": ra.uptime(), "reps": []}
    for r in range(7):
        out["reps"].append({"k4": pair(f"a2_4k{r}"), "thp": pair(f"a2_thp{r}", ex.TUN)})
        print(r, {c: {a: (round(out["reps"][-1][c][a]["user_s"], 2), round(out["reps"][-1][c][a]["sys_s"], 2)) for a in ("rho", "ic")} for c in ("k4", "thp")}, flush=True)
    out["load_after"] = ra.uptime()
    json.dump(out, open(os.path.join(HERE, "explore_thp_a2.json"), "w"), indent=2)
