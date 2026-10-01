"""A measurement session: one host, one reservation, interleaved arms.

The session reuses tools/isolated_bench.py for everything it already does
(AGENTS.md section 10 makes that tool mandatory for wall-clock numbers):
the exclusive lock, the CPU-mask checks including SMT siblings, the
settle-then-refuse preflight, and the eviction and restoration of other
threads.  ICMS adds the host capsule, per-run schedstat/steal/foreign-task
sampling, the graded isolation gate, interleaving, and structured records.

Plan.  Every (spec, workload) pair is an arm.  Each arm gets ``warmup``
unrecorded-for-comparison executions first (kept in the records with
``warmup: true``), then ``repetitions`` rounds.  In each round every arm runs
once; ``abab`` keeps the listed order, ``seeded_blocks`` shuffles each round
with a seed derived from the session id, so drift cannot line up with one
arm.  Listing the same spec twice gives an A/A pair.

Output (a new directory; an existing one is refused, never overwritten):

    session.json        plan, reservation, preflight, capsule hash, spec ids
    capsule.json        the host capsule
    records.jsonl       one icms.record/v1 per execution, in execution order
    specs/              a byte copy of every spec file used, with sha256 in session.json
    exec/<n>/           stdout, stderr and any producer files of execution n
"""
from __future__ import annotations

import contextlib
import importlib.util
import json
import os
import random
import shutil
import time
from typing import Any

from . import RECORD_SCHEMA, SESSION_SCHEMA, STANDARD, adapters
from .canonical import record_id, sha256_file, sha256_hex
from .environment import capsule as make_capsule, git_state, volatile_snapshot
from .execute import build_env, run_measured
from .gates import evaluate, level_at_least
from .registry import Registry
from .spec import REPO, load as load_spec


def _isolated_bench():
    path = os.path.join(REPO, "tools", "isolated_bench.py")
    spec = importlib.util.spec_from_file_location("isolated_bench", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class SessionError(RuntimeError):
    pass


def auto_cpus() -> set[int]:
    """The highest-numbered core this process may use, with all its SMT
    siblings: never CPU 0's core, never every CPU (there must be somewhere to
    move other work)."""
    ib = _isolated_bench()
    allowed = os.sched_getaffinity(0)
    for cpu in sorted(allowed, reverse=True):
        core = ib.smt_siblings(cpu) | {cpu}
        if 0 in core or not core <= allowed or core == allowed:
            continue
        return core
    raise SessionError(f"no core outside CPU 0's can be reserved from affinity {sorted(allowed)}")


def plan(arms: list[dict[str, Any]], session_seed: int) -> list[tuple[int, int, bool]]:
    """Return [(arm_index, round, is_warmup)] in execution order."""
    order: list[tuple[int, int, bool]] = []
    for i, arm in enumerate(arms):
        for w in range(arm["spec"]["measurement"]["warmup"]):
            order.append((i, -1 - w, True))
    rounds = max(a["spec"]["measurement"]["repetitions"] for a in arms)
    shuffle = any(a["spec"]["measurement"]["interleave"] == "seeded_blocks" for a in arms)
    rng = random.Random(session_seed)
    for r in range(rounds):
        block = [i for i, a in enumerate(arms) if r < a["spec"]["measurement"]["repetitions"]]
        if shuffle:
            rng.shuffle(block)
        order.extend((i, r, False) for i in block)
    return order


def run_session(spec_paths: list[str], out_dir: str, cpus: set[int], *, lock: str = "/tmp/crypto-bench.lock",
                wait: bool = False, settle: float = 2.0, max_other_cpu: float = 0.10, max_psi: float = 5.0,
                allow_busy: bool = False, session_label: str = "", log=print) -> dict[str, Any]:
    if os.path.exists(out_dir):
        raise SessionError(f"{out_dir} exists; a session never overwrites evidence")
    reg = Registry.load()
    loaded = [load_spec(p, reg) for p in spec_paths]
    problems = []
    arms = []
    for sp_path, sp in zip(spec_paths, loaded):
        ad = adapters.get(sp["spec"]["execution"]["adapter"])
        probs = ad.check(sp["spec"])
        problems += [f"{sp_path}: {x}" for x in probs]
        need = sp["spec"]["execution"]["cpus"]
        if need > len(cpus):
            problems.append(f"{sp_path}: needs {need} pinned CPUs, session reserved {len(cpus)}")
        for wl in sp["workloads"]:
            arms.append({"spec_path": sp_path, "spec": sp["spec"], "spec_id": sp["spec_id"], "workload": wl,
                         "adapter": ad})
    if problems:
        raise SessionError("refused before any run:\n  " + "\n  ".join(problems))

    started = time.time()
    session_id = "ICSESS1h" + sha256_hex({"specs": [a["spec_id"] for a in arms],
                                          "workloads": [a["workload"]["workload_id"] for a in arms],
                                          "started": repr(started), "label": session_label}, strict=False)[:12]
    os.makedirs(os.path.join(out_dir, "specs"))
    os.makedirs(os.path.join(out_dir, "exec"))
    spec_files = {}
    for p in dict.fromkeys(spec_paths):
        dst = os.path.join(out_dir, "specs", os.path.basename(p))
        shutil.copyfile(p, dst)
        spec_files[os.path.basename(p)] = sha256_file(dst)

    log(f"capturing host capsule")
    cap = make_capsule()
    with open(os.path.join(out_dir, "capsule.json"), "w") as fh:
        json.dump(cap, fh, indent=1, sort_keys=True)
    capsule_sha = sha256_file(os.path.join(out_dir, "capsule.json"))
    implementation = {"repo": git_state(REPO), "registry_sha256": reg.sha256,
                      "isolated_bench_sha256": sha256_file(os.path.join(REPO, "tools", "isolated_bench.py")),
                      "isolation_mode": "reserve-equivalent: the isolated_bench lock, preflight, CPU checks and "
                                        "eviction, imported from tools/isolated_bench.py; ICMS pins its own children",
                      "icms_sources": {f: sha256_file(os.path.join(os.path.dirname(__file__), f))
                                       for f in sorted(os.listdir(os.path.dirname(__file__))) if f.endswith(".py")}}

    ib = _isolated_bench()
    try:
        inherited = ib.check_cpus(cpus)
    except SystemExit as exc:
        shutil.rmtree(out_dir)
        raise SessionError(f"reservation refused: {exc}") from exc
    order = plan(arms, int(session_id[-8:], 16))
    session: dict[str, Any] = {
        "schema": SESSION_SCHEMA, "standard": STANDARD, "session_id": session_id, "label": session_label,
        "started_unix": started, "capsule_sha256": capsule_sha, "env_class_id": cap["env_class_id"],
        "implementation": implementation, "spec_files": spec_files,
        "arms": [{"index": i, "spec_path": a["spec_path"], "spec_id": a["spec_id"], "label": a["spec"].get("label"),
                  "role": a["spec"].get("role"), "workload_id": a["workload"]["workload_id"],
                  "adapter": a["adapter"].name} for i, a in enumerate(arms)],
        "plan": [{"arm": i, "round": r, "warmup": w} for i, r, w in order],
        "affinity_widened_from": inherited,
    }
    records_path = os.path.join(out_dir, "records.jsonl")
    lock_ctx = ib.locked(lock, wait)
    with lock_ctx:
        mine = {os.getpid()}
        args = type("A", (), {"settle": settle, "max_other_cpu": max_other_cpu, "max_psi": max_psi})()
        quiet, pre = True, None
        try:
            pre = ib.preflight(args, mine)
        except SystemExit as exc:
            quiet = False
            if not allow_busy:
                raise SessionError(f"preflight refused: {exc}") from exc
            pre = {"refused": str(exc)}
        worst = None
        if pre and "conditions" in pre:
            vals = [p["some"]["avg10"] for p in (pre["conditions"].get("psi_cpu"), pre["conditions"].get("psi_memory")) if p]
            worst = max(vals) if vals else None
        session["preflight"] = {"quiet": quiet, "detail": pre,
                                "other_cpu_seconds": (pre or {}).get("settle", {}).get("other_cpu_seconds"),
                                "psi_some_avg10_max": worst,
                                "limits": {"settle_s": settle, "max_other_cpu": max_other_cpu, "max_psi": max_psi}}
        moved = ib.evict(cpus, mine)
        session["reservation"] = {"cpus": sorted(cpus), "evicted": True, "threads_moved": len(moved),
                                  "left_on_reserved": ib.unmovable_on(cpus, mine)}
        session["volatile_start"] = volatile_snapshot()
        try:
            with open(records_path, "w") as rec_fh:
                for n, (ai, rnd, warm) in enumerate(order):
                    arm = arms[ai]
                    rec = _execute_one(n, ai, arm, rnd, warm, out_dir, cpus, cap, session, implementation, reg)
                    rec_fh.write(json.dumps(rec, sort_keys=True) + "\n")
                    rec_fh.flush()
                    log(f"[{n + 1}/{len(order)}] arm {ai} {'warmup' if warm else f'round {rnd}'}: "
                        f"{rec['outcome']['status']} {rec['execution']['wall_ns'] / 1e6:.1f} ms {rec['isolation']['earned_level']}"
                        + (f" (blocked: {','.join(rec['isolation']['blocking_next_level'])})" if rec['isolation']['blocking_next_level'] else ""))
        finally:
            ib.restore(moved)
    session["volatile_end"] = volatile_snapshot()
    session["finished_unix"] = time.time()
    session["records_sha256"] = sha256_file(records_path)
    with open(os.path.join(out_dir, "session.json"), "w") as fh:
        json.dump(session, fh, indent=1, sort_keys=True)
    return session


def _execute_one(n: int, arm_index: int, arm: dict[str, Any], rnd: int, warm: bool, out_dir: str, cpus: set[int],
                 cap: dict[str, Any], session: dict[str, Any], implementation: dict[str, Any],
                 reg: Registry) -> dict[str, Any]:
    spec, ad, wl = arm["spec"], arm["adapter"], arm["workload"]
    exec_dir = os.path.join(out_dir, "exec", f"{n:04d}")
    os.makedirs(exec_dir)
    ctx = adapters.Context(repo_root=REPO, exec_dir=exec_dir, out_dir=out_dir)
    cmd = ad.command(spec, wl, ctx)
    env, env_record = build_env(cmd["env"])
    want = spec["execution"]["cpus"]
    run_cpus = set(sorted(cpus)[:want])
    execution = run_measured(cmd["argv"], run_cpus, env, cmd["cwd"],
                             os.path.join(exec_dir, "stdout"), os.path.join(exec_dir, "stderr"),
                             timeout_s=spec["measurement"].get("timeout_seconds"),
                             memory_limit_mib=spec["measurement"].get("memory_limit_mib"),
                             interval=spec["measurement"]["sample_interval_ms"] / 1000.0)
    execution["env"] = env_record
    execution["argv"] = [os.path.relpath(a, REPO) if isinstance(a, str) and a.startswith(REPO + os.sep) else a
                         for a in execution["argv"]]
    execution["cwd"] = os.path.relpath(execution["cwd"], REPO) if execution.get("cwd") else None
    parsed = ad.parse(spec, wl, ctx, os.path.join(exec_dir, "stdout")) if not execution["exit"]["timed_out"] else \
        {"outcome": {"status": "timeout", "verified": False}}
    parsed = adapters.normalise(parsed, spec, reg)
    if execution["exit"]["returncode"] not in (0, None) and parsed.get("outcome", {}).get("status") == "complete":
        parsed["outcome"]["status"] = "error"
        parsed["outcome"]["reason"] = f"exit code {execution['exit']['returncode']}"
    gate = evaluate(execution, session, cap["stable"], spec["execution"]["threads"],
                    spec["measurement"].get("thresholds"))
    required = spec["measurement"]["isolation_required"]
    rec: dict[str, Any] = {
        "schema": RECORD_SCHEMA,
        "session_id": session["session_id"],
        "execution_index": n,
        "arm": arm_index,
        "round": rnd,
        "warmup": warm,
        "spec_id": arm["spec_id"],
        "spec_path": arm["spec_path"],
        "label": spec.get("label"),
        "role": spec.get("role"),
        "workload_id": wl["workload_id"],
        "workload": wl["record"],
        "adapter": ad.name,
        "window": spec["measurement"]["window"],
        "unit": spec["accounting"]["unit"],
        "implementation": {"repo_commit": (implementation.get("repo") or {}).get("commit"),
                           "repo_dirty": (implementation.get("repo") or {}).get("dirty"),
                           "producer": ad.provenance(spec, ctx)},
        "environment": {"env_class_id": cap["env_class_id"], "capsule_sha256": session["capsule_sha256"]},
        "execution": execution,
        "isolation": {**gate, "required": required, "wall_admissible": level_at_least(gate["earned_level"], required)},
        "outcome": parsed.get("outcome"),
        "units": parsed.get("units"),
        "reference": parsed.get("reference"),
        "metrics": parsed.get("metrics"),
        "phases": parsed.get("phases"),
        "windows": parsed.get("windows"),
        "consistency": parsed.get("consistency") or [],
        "producer": parsed.get("producer"),
    }
    if rec.get("windows") and "whole_process" in rec["windows"]:
        rec["windows"]["whole_process"]["wall_ns"] = execution["wall_ns"]
    rec["record_id"] = record_id(rec)
    return rec
