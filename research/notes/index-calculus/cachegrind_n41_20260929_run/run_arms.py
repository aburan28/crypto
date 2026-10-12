#!/usr/bin/env python3
"""Measurement stages for RESEARCH_CACHEGRIND_N41_20260929.md (frozen cell n=41, L=1,024, K=255).

usage: run_arms.py <stage>          (env MEM_CALIB=<path to mem_calib binary>)
  native        5 interleaved reps of rho R3 and IC, uninstrumented, pinned to core 0
  cachegrind    both arms x LL in {2, 8, 32} MiB under valgrind --tool=cachegrind (--branch-sim=yes)
  controls      simulator sanity controls (mem_calib under cachegrind) and their native timings
  interference  3 reps of each arm alone vs with three memory antagonists on cores 1-3

All outputs land in work/ next to this script; each stage writes <stage>.json.
"""
import concurrent.futures
import json
import os
import resource
import shutil
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
EX = os.path.join(ROOT, "target", "release", "examples")
RHO = os.path.join(EX, "koblitz_rho_batch_ks_strong")
IC = os.path.join(EX, "koblitz_orbit_dlp_fast")
MEM = os.environ.get("MEM_CALIB", "")
SWEEP_CELL = os.path.join(HERE, "..", "strong_rho_sweep_20260929_run", "cell_n41_L1024_K255")
WORK = os.path.join(HERE, "work")
N, L, K = 41, 1024, 255
CORPUS, BATCH_SEED = "n41-strong-sweep-L1024-v1", "531310"
CG_COMMON = ["--tool=cachegrind", "--cache-sim=yes", "--branch-sim=yes",
             "--I1=32768,8,64", "--D1=32768,8,64"]
MIB = 1024 * 1024
LL_SIZES = (2 * MIB, 8 * MIB, 32 * MIB)


def uptime():
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def rho_env():
    env = dict(os.environ)
    env.update(KIC_RHO_RUNG="3", KIC_RHO_LANES="32", KIC_RHO_DP_BITS="4", KIC_RHO_BATCH_CORPUS=CORPUS)
    return env


def rho_cmd():
    return [RHO, str(N), "0", "signed_frobenius", str(L), BATCH_SEED]


def ic_cmd(out_jsonl):
    return [IC, f"construct:{N}:0:{K}", os.path.join(WORK, "scalars.txt"), "7", out_jsonl]


def parse_rho(path):
    fixtures, summary = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                obj = json.loads(line)
                if obj.get("kind") == "rho_ks_batch_fixture":
                    fixtures.append(obj)
                elif obj.get("kind") == "rho_ks_batch_summary":
                    summary = obj
    return fixtures, summary


def scalars():
    return [int(x) for x in open(os.path.join(WORK, "scalars.txt")).read().split()]


def rho_ok(path, rc):
    fixtures, summary = parse_rho(path)
    return (rc == 0 and len(fixtures) == L and summary is not None and summary["all_verified"] is True
            and all(f["verified"] and f["recovered_fixture_scalar"] == f["published_fixture_scalar"] for f in fixtures)
            and [int(f["published_fixture_scalar"]) for f in fixtures] == scalars()), \
        (summary or {}).get("total_walk_steps")


def ic_ok(path, rc):
    try:
        s = json.load(open(path))
    except Exception:  # noqa: BLE001
        return False, None
    return rc == 0 and s.get("targets_solved") == L and s.get("targets_failed") == 0, s


def run_timed(cmd, out, err, env=None, pin=None):
    """Run cmd; return rc, wall_s, cpu_s (user+sys of the child), maxrss_bytes."""
    if pin is not None:
        cmd = ["taskset", "-c", str(pin)] + cmd
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    with open(out, "w") as fo, open(err, "w") as fe:
        rc = subprocess.run(cmd, stdout=fo, stderr=fe, env=env).returncode
    wall = time.perf_counter() - t0
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)
    return rc, wall, cpu, r1.ru_maxrss * 1024


def prepare():
    os.makedirs(WORK, exist_ok=True)
    dst = os.path.join(WORK, "scalars.txt")
    if not os.path.exists(dst):
        shutil.copy(os.path.join(SWEEP_CELL, "scalars.txt"), dst)


def one_pair(tag, pin, env_extra_cmd_prefix=None):
    """Run rho then IC once; return a record."""
    rec = {}
    rc, wall, cpu, rss = run_timed(rho_cmd(), f"{WORK}/rho_{tag}.jsonl", f"{WORK}/rho_{tag}.stderr", rho_env(), pin)
    ok, steps = rho_ok(f"{WORK}/rho_{tag}.jsonl", rc)
    rec["rho"] = {"rc": rc, "wall_s": wall, "cpu_s": cpu, "gates_G1": ok, "total_walk_steps": steps}
    rc, wall, cpu, rss = run_timed(ic_cmd(f"{WORK}/ic_{tag}.jsonl"), f"{WORK}/ic_{tag}.summary.json",
                                   f"{WORK}/ic_{tag}.stderr", None, pin)
    ok, s = ic_ok(f"{WORK}/ic_{tag}.summary.json", rc)
    rec["ic"] = {"rc": rc, "wall_s": wall, "cpu_s": cpu, "gates_G1": ok,
                 "peak_rss_bytes": (s or {}).get("peak_rss_bytes") if isinstance(s, dict) else None}
    return rec


def stage_native(reps=5):
    prepare()
    out = {"stage": "native", "reps": [], "load_before": uptime()}
    for r in range(reps):
        out["reps"].append(one_pair(f"native{r}", pin=0))
    out["load_after"] = uptime()
    json.dump(out, open(os.path.join(HERE, "native.json"), "w"), indent=2)
    for arm in ("rho", "ic"):
        c = [x[arm]["cpu_s"] for x in out["reps"]]
        print(arm, "cpu_s", [round(v, 3) for v in c], "spread", round(max(c) / min(c) - 1, 4),
              "G1", all(x[arm]["gates_G1"] for x in out["reps"]))


def parse_cg(path):
    events, summary = None, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("events:"):
                events = line.split()[1:]
            elif line.startswith("summary:"):
                summary = [int(x) for x in line.split()[1:]]
    return dict(zip(events, summary))


def cg_job(arm, ll):
    tag = f"cg_{arm}_ll{ll // MIB}M"
    cmd = ["valgrind"] + CG_COMMON + [f"--LL={ll},16,64", f"--cachegrind-out-file={WORK}/{tag}.out"]
    if arm == "rho":
        cmd += rho_cmd()
        env, out = rho_env(), f"{WORK}/{tag}.stdout.jsonl"
    else:
        out = f"{WORK}/{tag}.summary.json"
        cmd += ic_cmd(f"{WORK}/{tag}.stdout.jsonl")
        env = None
    rc, wall, _, _ = run_timed(cmd, out, f"{WORK}/{tag}.stderr", env)
    ok = rho_ok(out, rc)[0] if arm == "rho" else ic_ok(out, rc)[0]
    return {"tag": tag, "arm": arm, "ll_bytes": ll, "rc": rc, "wall_s": wall, "gates_G1": ok,
            "counts": parse_cg(f"{WORK}/{tag}.out")}


def stage_cachegrind():
    prepare()
    jobs = [(arm, ll) for ll in LL_SIZES for arm in ("ic", "rho")]
    out = {"stage": "cachegrind", "cg_args": CG_COMMON + ["--LL=<size>,16,64"], "load_before": uptime(), "runs": []}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for res in ex.map(lambda a: cg_job(*a), jobs):
            out["runs"].append(res)
            print(res["tag"], res["gates_G1"], res["counts"]["Ir"])
    out["load_after"] = uptime()
    json.dump(out, open(os.path.join(HERE, "cachegrind.json"), "w"), indent=2)


# Simulator sanity controls: (label, mode, W bytes, k)
CONTROLS = [("chase_512K_k1", "chase", 512 * 1024, 1), ("chase_8M_k1", "chase", 8 * MIB, 1),
            ("chase_256M_k1", "chase", 256 * MIB, 1), ("chase_256M_k10", "chase", 256 * MIB, 10),
            ("seq_256M", "seq", 256 * MIB, 0)]
CTL_STEPS = 4_000_000  # group steps (chase) or passes*lines handled below
CTL_S0 = 1000


def ctl_args(mode, w, k, steps):
    if mode == "chase":
        return [MEM, "chase", str(w), str(k), str(steps), "1", "0"]
    return [MEM, "seq", str(w), str(steps), "1"]  # for seq, "steps" is the number of passes


def cg_ctl(label, mode, w, k, steps, ll):
    tag = f"ctl_{label}_s{steps}_ll{ll // MIB}M"
    cmd = ["valgrind"] + CG_COMMON + [f"--LL={ll},16,64", f"--cachegrind-out-file={WORK}/{tag}.out"] + ctl_args(mode, w, k, steps)
    rc, _, _, _ = run_timed(cmd, f"{WORK}/{tag}.stdout", f"{WORK}/{tag}.stderr")
    return {"label": label, "steps": steps, "ll_bytes": ll, "rc": rc, "counts": parse_cg(f"{WORK}/{tag}.out")}


def stage_controls():
    prepare()
    out = {"stage": "controls", "load_before": uptime(), "steps": CTL_STEPS, "s0": CTL_S0, "cachegrind": [], "native": []}
    sv = lambda mode: (CTL_STEPS, CTL_S0) if mode == "chase" else (2, 1)
    jobs = [(lab, mode, w, k, s, ll) for (lab, mode, w, k) in CONTROLS for ll in LL_SIZES[:2] for s in sv(mode)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for res in ex.map(lambda a: cg_ctl(*a), jobs):
            out["cachegrind"].append(res)
            print(res["label"], res["steps"], res["ll_bytes"] // MIB, res["counts"]["Ir"])
    # native timings, sequential, pinned; baseline = same k at 16 KiB
    for lab, mode, w, k in CONTROLS + [("base_k1", "chase", 16384, 1), ("base_k10", "chase", 16384, 10)]:
        if mode == "chase":
            cmd = [MEM, "chase", str(w), str(k), str(CTL_STEPS), "5", str(CTL_STEPS // 10)]
        else:
            cmd = [MEM, "seq", str(w), "2", "5"]
        res = subprocess.run(["taskset", "-c", "0"] + cmd, capture_output=True, text=True, check=True).stdout
        out["native"].append({"label": lab, **json.loads(res.strip().splitlines()[-1])})
        print("native", lab, out["native"][-1])
    out["load_after"] = uptime()
    json.dump(out, open(os.path.join(HERE, "controls.json"), "w"), indent=2)


def stage_interference(reps=3):
    prepare()
    out = {"stage": "interference", "reps": [], "load_before": uptime()}
    for r in range(reps):
        alone = one_pair(f"alone{r}", pin=0)
        ants = [subprocess.Popen(["taskset", "-c", str(c), MEM, "antagonist", str(256 * MIB), "8", "120"],
                                 stdout=subprocess.DEVNULL) for c in (1, 2, 3)]
        time.sleep(8)  # let each antagonist build its 256 MiB cycle set
        loaded = one_pair(f"loaded{r}", pin=0)
        for a in ants:
            a.terminate()
        for a in ants:
            a.wait()
        out["reps"].append({"alone": alone, "loaded": loaded})
        print(r, {a: (round(alone[a]["cpu_s"], 2), round(loaded[a]["cpu_s"], 2)) for a in ("rho", "ic")})
    out["load_after"] = uptime()
    json.dump(out, open(os.path.join(HERE, "interference.json"), "w"), indent=2)


if __name__ == "__main__":
    if not MEM or not os.path.exists(MEM):
        sys.exit("set MEM_CALIB to the mem_calib binary")
    {"native": stage_native, "cachegrind": stage_cachegrind, "controls": stage_controls,
     "interference": stage_interference}[sys.argv[1]]()
