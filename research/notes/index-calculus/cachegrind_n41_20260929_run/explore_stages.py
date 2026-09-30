#!/usr/bin/env python3
"""Unregistered follow-up stages of Amendment A1 in RESEARCH_CACHEGRIND_N41_20260929.md.

usage: explore_stages.py <stage>        (env MEM_CALIB=<mem_calib binary>)
  replicate     E3: 10 more interleaved native repetitions of (rho, IC), core 0, as in the registered stage
  thp           E2: 7 interleaved repetitions of {rho, IC} x {4 KiB pages, GLIBC_TUNABLES=glibc.malloc.hugetlb=1};
                    also samples AnonHugePages from /proc/<pid>/smaps_rollup; plus a 256 MiB pointer-chase
                    control with and without the tunable
  interference  E1: positive/negative control for the antagonist test (pointer chase over 256 MiB and over
                    16 KiB, alone vs with the same three antagonists)
Outputs: explore_<stage>.json next to this script; per-run logs in work/.
"""
import importlib.util
import json
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ra", os.path.join(HERE, "run_arms.py"))
ra = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ra)
MEM, WORK, MIB = ra.MEM, ra.WORK, ra.MIB
TUN = {"GLIBC_TUNABLES": "glibc.malloc.hugetlb=1"}


def anon_huge_kb(pid):
    try:
        for line in open(f"/proc/{pid}/smaps_rollup"):
            if line.startswith("AnonHugePages:"):
                return int(line.split()[1])
    except OSError:
        pass
    return None


def run_polled(cmd, out, err, env, pin=0):
    """Run cmd pinned; return rc, wall, cpu, max AnonHugePages (kB) seen while it ran."""
    full = ["taskset", "-c", str(pin)] + cmd
    t0 = time.perf_counter()
    with open(out, "w") as fo, open(err, "w") as fe:
        p = subprocess.Popen(full, stdout=fo, stderr=fe, env=env)
        best = 0
        while True:
            done, status, ru = os.wait4(p.pid, os.WNOHANG)
            if done:
                break
            v = anon_huge_kb(p.pid)
            if v:
                best = max(best, v)
            time.sleep(0.25)
    wall = time.perf_counter() - t0
    rc = os.waitstatus_to_exitcode(status)
    return rc, wall, ru.ru_utime + ru.ru_stime, best


def pair(tag, env_extra=None, pin=0):
    ra.prepare()
    rec = {}
    renv = ra.rho_env()
    if env_extra:
        renv.update(env_extra)
    rc, wall, cpu, huge = run_polled(ra.rho_cmd(), f"{WORK}/rho_{tag}.jsonl", f"{WORK}/rho_{tag}.stderr", renv, pin)
    ok, steps = ra.rho_ok(f"{WORK}/rho_{tag}.jsonl", rc)
    rec["rho"] = {"rc": rc, "wall_s": wall, "cpu_s": cpu, "gates_G1": ok, "max_anon_huge_kb": huge}
    ienv = dict(os.environ)
    if env_extra:
        ienv.update(env_extra)
    rc, wall, cpu, huge = run_polled(ra.ic_cmd(f"{WORK}/ic_{tag}.jsonl"), f"{WORK}/ic_{tag}.summary.json",
                                     f"{WORK}/ic_{tag}.stderr", ienv, pin)
    ok, s = ra.ic_ok(f"{WORK}/ic_{tag}.summary.json", rc)
    rec["ic"] = {"rc": rc, "wall_s": wall, "cpu_s": cpu, "gates_G1": ok, "max_anon_huge_kb": huge,
                 "peak_rss_bytes": s.get("peak_rss_bytes") if isinstance(s, dict) else None}
    return rec


def chase(w, k, steps, reps, env_extra=None, pin=0):
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    out = subprocess.run(["taskset", "-c", str(pin), MEM, "chase", str(w), str(k), str(steps), str(reps), str(steps // 10)],
                         capture_output=True, text=True, check=True, env=env).stdout
    return json.loads(out.strip().splitlines()[-1])


def stage_replicate(reps=10):
    out = {"stage": "replicate", "load_before": ra.uptime(), "reps": [pair(f"rep{r}") for r in range(reps)]}
    out["load_after"] = ra.uptime()
    json.dump(out, open(os.path.join(HERE, "explore_replicate.json"), "w"), indent=2)
    for a in ("rho", "ic"):
        c = [x[a]["cpu_s"] for x in out["reps"]]
        print(a, [round(v, 3) for v in c], "spread", round(max(c) / min(c) - 1, 3))


def stage_thp(reps=7):
    out = {"stage": "thp", "load_before": ra.uptime(), "reps": [], "chase_control": []}
    for r in range(reps):
        out["reps"].append({"k4": pair(f"thp4k{r}"), "thp": pair(f"thpon{r}", TUN)})
        print(r, {c: {a: round(out["reps"][-1][c][a]["cpu_s"], 2) for a in ("rho", "ic")} for c in ("k4", "thp")},
              "AnonHuge kB", {a: out["reps"][-1]["thp"][a]["max_anon_huge_kb"] for a in ("rho", "ic")})
    for r in range(3):
        out["chase_control"].append({"k4": chase(256 * MIB, 1, 4_000_000, 3),
                                     "thp": chase(256 * MIB, 1, 4_000_000, 3, TUN)})
    out["load_after"] = ra.uptime()
    json.dump(out, open(os.path.join(HERE, "explore_thp.json"), "w"), indent=2)


def stage_interference():
    out = {"stage": "interference_controls", "load_before": ra.uptime(), "reps": []}
    for r in range(3):
        rec = {"alone": {"pos": chase(256 * MIB, 1, 6_000_000, 3), "neg": chase(16384, 1, 200_000_000, 3)}}
        ants = [subprocess.Popen(["taskset", "-c", str(c), MEM, "antagonist", str(256 * MIB), "8", "120"],
                                 stdout=subprocess.DEVNULL) for c in (1, 2, 3)]
        time.sleep(8)
        rec["loaded"] = {"pos": chase(256 * MIB, 1, 6_000_000, 3), "neg": chase(16384, 1, 200_000_000, 3)}
        for a in ants:
            a.terminate()
        for a in ants:
            a.wait()
        out["reps"].append(rec)
        print(r, {k: {t: round(rec[k][t]["ns_per_step_median"], 2) for t in ("pos", "neg")} for k in rec})
    out["load_after"] = ra.uptime()
    json.dump(out, open(os.path.join(HERE, "explore_interference.json"), "w"), indent=2)


if __name__ == "__main__":
    if not MEM or not os.path.exists(MEM):
        sys.exit("set MEM_CALIB")
    {"replicate": stage_replicate, "thp": stage_thp, "interference": stage_interference}[sys.argv[1]]()
