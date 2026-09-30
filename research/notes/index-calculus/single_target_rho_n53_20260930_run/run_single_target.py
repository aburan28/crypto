#!/usr/bin/env python3
"""Measurement stages for RESEARCH_SINGLE_TARGET_STRONG_RHO_N53_20260930.md.

usage: run_single_target.py <stage>
  timing     9 randomized blocks of D (direct, 4 threads), P (panel rho), S (strong rho), whole process
  seeds      P and S under 16 walk seeds each, one run per seed, single thread
  callgrind  one callgrind run each of D, P, S (whole-process Ir)
Outputs: <stage>.json plus raw logs in work/.  Env: none (paths derived from the repository).
"""
import hashlib
import json
import os
import pathlib
import random
import resource
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import run_koblitz_stage106_shard_routing as s106  # noqa: E402  (reuses the panel's environment and parsers)
import run_koblitz_stage42_n53_same_target as s42  # noqa: E402

EX = os.path.join(ROOT, "target", "release", "examples")
DIRECT = os.path.join(EX, "koblitz_rank_fixture")
PANEL_RHO = os.path.join(EX, "koblitz_rho_fixture")
STRONG = os.path.join(EX, "koblitz_rho_batch_ks_strong")
WORK = os.path.join(HERE, "work")
SCALAR = s42.EXPLICIT_SCALAR            # 476811900269
PANEL_SEED = s42.RHO_BATCH_SEED         # 16969290756857938023
STRONG_SEED0 = 531310
BLOCKS, NSEEDS, RNG_SEED = 9, 16, 20260930


def seed_for(arm, i):
    if i == 0:
        return PANEL_SEED if arm == "P" else STRONG_SEED0
    digest = hashlib.sha256(f"SINGLE-TARGET-RHO-N53|{arm}|{i}".encode()).digest()
    return int.from_bytes(digest[:8], "little")


def cmd_D():
    return [DIRECT, "53", "0", "1", "128", str(s42.DIRECT_SEED), "signed_expanded", "independent",
            "pair_pair_parallel_4096", "1", str(SCALAR)]


def cmd_P(seed):
    return [PANEL_RHO, "53", "0", "signed_frobenius", "1", "packed", str(seed), str(SCALAR)]


def cmd_S(seed):
    return [STRONG, "53", "0", "signed_frobenius", "1", str(seed)]


def env_for(arm):
    if arm == "D":
        return s106.common_environment()
    env = s106.custody.safe_child_environment()
    if arm == "S":
        env.update(KIC_RHO_RUNG="3", KIC_RHO_LANES="32", KIC_RHO_DP_BITS="4", KIC_RHO_EXPLICIT_SCALAR=str(SCALAR))
    return env


def pinned(arm, cmd):
    return ["taskset", "-c", "0-3" if arm == "D" else "0"] + cmd


def run_one(arm, cmd, tag, wrap=None):
    os.makedirs(WORK, exist_ok=True)
    out, err = f"{WORK}/{tag}.stdout", f"{WORK}/{tag}.stderr"
    full = (wrap or []) + pinned(arm, cmd)
    t0 = time.perf_counter()
    with open(out, "w") as fo, open(err, "w") as fe:
        p = subprocess.Popen(full, stdout=fo, stderr=fe, env=env_for(arm))
        _, status, ru = os.wait4(p.pid, 0)
    wall = time.perf_counter() - t0
    return {"tag": tag, "arm": arm, "rc": os.waitstatus_to_exitcode(status), "wall_s": wall,
            "user_s": ru.ru_utime, "sys_s": ru.ru_stime, "cpu_s": ru.ru_utime + ru.ru_stime,
            "maxrss_kb": ru.ru_maxrss, "stdout": out, "stderr": err}


def jsonl(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if line.startswith("{"):
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return rows


def check(rec):
    """Completion criteria; returns (ok, extra fields)."""
    rows = jsonl(rec["stdout"])
    if rec["rc"] != 0:
        return False, {}
    if rec["arm"] == "D":
        try:
            obs = s106.stage61.direct_observation(pathlib.Path(rec["stdout"]), "FULL_RANK")
        except Exception as exc:  # noqa: BLE001 - recorded as a failed gate
            return False, {"error": str(exc)}
        return True, {"admitted_relations": obs["summary"].get("admitted_relations")}
    if rec["arm"] == "P":
        row = next((r for r in rows if r.get("kind") == "rho_public_fixture"), None)
        ok = bool(row and row.get("verified") and row.get("recovered_fixture_scalar") == SCALAR)
        return ok, {"walk_steps": row and row["walk_steps"], "ideal_steps": row and row["ideal_steps"],
                    "table_entries": row and row["table_entries"]}
    fx = [r for r in rows if r.get("kind") == "rho_ks_batch_fixture"]
    ok = len(fx) == 1 and fx[0].get("verified") and fx[0].get("recovered_fixture_scalar") == SCALAR \
        and fx[0].get("published_fixture_scalar") == SCALAR
    return bool(ok), {"walk_steps": fx and fx[0]["walk_steps"], "ideal_steps": fx and fx[0]["ideal_independent_steps"]}


def stage_timing():
    rng = random.Random(RNG_SEED)
    out = {"stage": "timing", "load_before": subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip(), "blocks": []}
    for b in range(BLOCKS):
        order = ["D", "P", "S"]
        rng.shuffle(order)
        blk = {"order": order}
        for arm in order:
            cmd = {"D": cmd_D(), "P": cmd_P(seed_for("P", 0)), "S": cmd_S(seed_for("S", 0))}[arm]
            rec = run_one(arm, cmd, f"t{b}_{arm}")
            rec["gate"], rec["extra"] = check(rec)
            blk[arm] = rec
        out["blocks"].append(blk)
        print(b, order, {a: (round(blk[a]["wall_s"], 2), round(blk[a]["cpu_s"], 2), blk[a]["gate"]) for a in "DPS"}, flush=True)
    out["load_after"] = subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()
    json.dump(out, open(os.path.join(HERE, "timing.json"), "w"), indent=2)


def stage_seeds():
    out = {"stage": "seeds", "runs": []}
    for arm, cmdf in (("P", cmd_P), ("S", cmd_S)):
        for i in range(NSEEDS):
            seed = seed_for(arm, i)
            rec = run_one(arm, cmdf(seed), f"s_{arm}{i}")
            rec["seed"], rec["index"] = seed, i
            rec["gate"], rec["extra"] = check(rec)
            out["runs"].append(rec)
            print(arm, i, round(rec["wall_s"], 3), rec["gate"], rec["extra"].get("walk_steps"), flush=True)
    json.dump(out, open(os.path.join(HERE, "seeds.json"), "w"), indent=2)


def stage_callgrind():
    out = {"stage": "callgrind", "runs": []}
    for arm, cmd in (("D", cmd_D()), ("P", cmd_P(seed_for("P", 0))), ("S", cmd_S(seed_for("S", 0)))):
        rec = run_one(arm, cmd, f"cg_{arm}", wrap=["valgrind", "--tool=callgrind", f"--callgrind-out-file={WORK}/cg_{arm}.out"])
        rec["gate"], rec["extra"] = check(rec)
        ir = None
        for line in open(rec["stderr"]):
            if "Collected :" in line:
                ir = int(line.split("Collected :")[1].split()[0])
        rec["Ir"] = ir
        out["runs"].append(rec)
        print(arm, rec["rc"], rec["gate"], ir, round(rec["wall_s"], 1), flush=True)
    json.dump(out, open(os.path.join(HERE, "callgrind.json"), "w"), indent=2)


if __name__ == "__main__":
    {"timing": stage_timing, "seeds": stage_seeds, "callgrind": stage_callgrind}[sys.argv[1]]()
