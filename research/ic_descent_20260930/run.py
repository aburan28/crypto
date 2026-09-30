#!/usr/bin/env python3
"""Ledger §22's runs, in the order PROTOCOL.md declares them.

    python3 run.py manifest   # host manifest -> host.json
    python3 run.py inputs     # check §20's parameter files against inputs.sha256
    python3 run.py aa         # v2: baseline against a byte-identical copy, M1, nine sizes
    python3 run.py main       # 36 files x 5 rounds, baseline then candidate
    python3 run.py control1   # ic workflow == ic price, candidate, M1 at every size
    python3 run.py rho        # batch rho re-priced on the baseline at n = 41, M1
    python3 run.py threads    # many threads, n = 41, M1, 3 ABAB rounds
    python3 run.py probe      # the descent probe on both libraries, six sizes
    python3 run.py all        # every step of the run's version, in its order

The binaries are named by the environment: IC_BASELINE and IC_CANDIDATE
for `ic`, PROBE_BASELINE and PROBE_CANDIDATE for
`examples/koblitz_descent_prices`, with IC_BASELINE_COMMIT and
IC_CANDIDATE_COMMIT naming the commits they were built from.  None is in
the tree; the manifest records each one's sha256.

**v1** (the first run, `runs/`): every single-thread process runs with
RAYON_NUM_THREADS=1 under `taskset -c 2`; the thread check with four
threads under `taskset -c 0-3`.

**v2** (`IC_ISOLATE=1 IC_RUNS=<dir>`, PROTOCOL.md v2): every timed
process runs through `tools/isolated_bench.py run --wait`, which locks
out other benchmarks, refuses a busy machine, reserves the CPUs and
records the conditions beside the report (`*.isolation.jsonl`).  One
thread on CPU 2; the thread check at three threads on CPUs 1-3, since
the tool will not reserve all four.  A refused start is logged to
`refusals.log` and retried after 15 s.  A pair with a contended (or
failed) process is kept and run again in the same slot, up to twice;
the analysis uses the first clean pair of each slot.  `aa` needs
IC_BASELINE_COPY, a byte-identical copy of the baseline `ic`.

Outputs go to IC_RUNS (default runs/), one file per process, and an
existing file is never overwritten: a rerun resumes where the last one
stopped.  Nothing is dropped; a failed run keeps its report and stderr.
A pin mismatch (counts or recovered logarithms differing between the
arms of a pair) stops the run, as declared.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
ISOLATE = os.environ.get("IC_ISOLATE") == "1"
TOOL = ROOT / "tools" / "isolated_bench.py"
S20 = ROOT / "research" / "ic_exponent_20260926"
# §20's nine sizes, in its order (by the size of r).
SIZES = [(1, 19), (1, 23), (1, 45), (0, 37), (1, 43), (1, 47), (0, 41), (0, 53), (0, 61)]
PROBED = [(1, 19), (1, 23), (1, 45), (0, 37), (1, 43), (0, 41)]
SETS = (1, 2, 3, 4)
ROUNDS = 5
SPREAD_LIMIT = 1.25
RETRIES = 2 if ISOLATE else 0
REFUSAL_WAIT_S, REFUSAL_LIMIT = 15, 60
ARMS = {"baseline": "IC_BASELINE", "candidate": "IC_CANDIDATE", "copy": "IC_BASELINE_COPY"}
PROBES = {"baseline": "PROBE_BASELINE", "candidate": "PROBE_CANDIDATE"}
# (environment, CPUs): pinned with taskset in v1, reserved by the tool in v2.
ONE = ({**os.environ, "RAYON_NUM_THREADS": "1"}, "2")
FOUR = ({**os.environ, "RAYON_NUM_THREADS": "4"}, "0-3")
THREE = ({**os.environ, "RAYON_NUM_THREADS": "3"}, "1-3")


def sh(cmd: list[str]) -> str:
    """A command's output, or "" when the tool is not installed."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout.strip()
    except FileNotFoundError:
        return ""


def binary(var: str) -> Path:
    path = os.environ.get(var)
    if not path or not Path(path).exists():
        raise SystemExit(f"set {var} to the binary")
    return Path(path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def params(a: int, n: int, j: int) -> Path:
    return S20 / "runs" / f"k{a}n{n}" / "measure" / f"M{j}.params.json"


def record_path(out: Path) -> Path:
    """Where the isolation record of the process that wrote `out` goes."""
    name = out.name
    for suffix in (".price.json", ".json"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return out.with_name(name + ".isolation.jsonl")


def clean(out: Path) -> bool:
    """v2: the process ran uncontended and exited 0.  v1 has no records."""
    if not ISOLATE:
        return True
    rec = record_path(out)
    if not rec.exists():
        return False
    run = json.loads(rec.read_text().splitlines()[-1])["run"]
    return run["exit_status"] == 0 and not run["contended"]


def launch(cmd: list[str], env, out: Path, stdout, stderr: Path) -> None:
    """Run one timed process: under taskset (v1) or through the tool (v2)."""
    environ, cpus = env
    if not ISOLATE:
        with open(stderr, "w") as e:
            subprocess.run(["taskset", "-c", cpus, *cmd], env=environ, stdout=stdout, stderr=e, check=False)
        return
    rec = record_path(out)
    label = str(out.relative_to(RUNS))
    for attempt in range(REFUSAL_LIMIT):
        tmp = stderr.with_suffix(".stderr.tmp")
        with open(tmp, "w") as e:
            subprocess.run([sys.executable, str(TOOL), "run", "--wait", "--cpus", cpus, "--out", str(rec),
                            "--label", label, "--", *cmd], env=environ, stdout=stdout, stderr=e, check=False)
        if rec.exists():
            tmp.replace(stderr)
            return
        # Refused before starting: the machine was busy.  Log it and wait.
        with open(RUNS / "refusals.log", "a") as log:
            log.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {label} attempt {attempt + 1}: "
                      f"{tmp.read_text().strip()}\n")
        tmp.unlink()
        if stdout not in (None, subprocess.DEVNULL) and hasattr(stdout, "seek"):
            stdout.seek(0)
            stdout.truncate()
        time.sleep(REFUSAL_WAIT_S)
    raise SystemExit(f"{label}: the machine stayed busy for {REFUSAL_LIMIT} tries; stopped (resumable)")


def manifest() -> None:
    out = (HERE / "host.json") if RUNS == (HERE / "runs").resolve() else (RUNS / "host.json")
    if out.exists():
        print(f"{out} exists; not overwritten")
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    cpu = {}
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if ":" in line:
            k, v = (s.strip() for s in line.split(":", 1))
            if k in ("model name", "flags") and k not in cpu:
                cpu[k] = v
    flags = cpu.get("flags", "").split()
    wanted = [f for f in flags if f in ("popcnt", "avx2", "pclmulqdq", "bmi2", "gfni", "vpclmulqdq")
              or f.startswith("avx512")]
    mem = next((l for l in Path("/proc/meminfo").read_text().splitlines() if l.startswith("MemTotal")), "")
    arms = [a for a in ARMS if os.environ.get(ARMS[a])]
    doc = {
        "commit": sh(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        "tree_status": sh(["git", "-C", str(ROOT), "status", "--porcelain"]),
        "rustc": sh(["rustc", "--version"]),
        "cpu_model": cpu.get("model name"),
        "cpu_flags_relevant": wanted,
        "logical_cores": os.cpu_count(),
        "memory": mem,
        "os": platform.platform(),
        "arch": platform.machine(),
        "binaries": {
            arm: {"ic_sha256": sha256(binary(ARMS[arm])),
                  **({"probe_sha256": sha256(binary(PROBES[arm]))} if arm in PROBES else {}),
                  "built_from": os.environ.get(f"{ARMS[arm]}_COMMIT")}
            for arm in arms
        },
        "build": "cargo build --release --bin ic --example koblitz_descent_prices, clean tree, rustc above",
        "uptime": sh(["uptime"]),
        "hardware_class": "one x86-64 cloud container; no claim for Arm64, GPUs or other hosts",
    }
    if ISOLATE:
        doc["pinning"] = ("tools/isolated_bench.py run --wait: RAYON_NUM_THREADS=1 on CPU 2; the thread check "
                          "RAYON_NUM_THREADS=3 on CPUs 1-3; the tool's defaults (settle 2 s, other processes "
                          "at most 0.10 CPUs, PSI some avg10 at most 5)")
        doc["isolation_tool"] = {"path": "tools/isolated_bench.py", "sha256": sha256(TOOL)}
        doc["smt_siblings"] = {c: Path(f"/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list")
                               .read_text().strip() for c in range(os.cpu_count() or 0)}
    else:
        doc["pinning"] = ("RAYON_NUM_THREADS=1 under taskset -c 2; the thread check RAYON_NUM_THREADS=4 "
                          "under taskset -c 0-3")
    out.write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))


def inputs() -> None:
    bad = []
    for line in (HERE / "inputs.sha256").read_text().splitlines():
        digest, rel = line.split()
        if sha256(S20 / "runs" / rel) != digest:
            bad.append(rel)
    if bad:
        raise SystemExit(f"inputs changed: {bad}")
    print("inputs: all 36 match inputs.sha256")


def price(arm: str, p: Path, out: Path, extra: list[str], env=ONE) -> dict:
    """Run `ic price` once into `out` unless it exists; return the report."""
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [str(binary(ARMS[arm])), "price", "--params", str(p), "--json", "--out", str(out), *extra]
        launch(cmd, env, out, subprocess.DEVNULL, out.with_suffix(".stderr"))
    return json.loads(out.read_text())


def pinned(a: dict, b: dict) -> bool:
    return (a.get("status") == b.get("status") == "complete"
            and a.get("counts") == b.get("counts") and a.get("recovered") == b.get("recovered")
            and a.get("all_verified") is True and b.get("all_verified") is True)


def rounds(d: Path, p: Path, tag: str, extra: list[str], count: int, env=ONE,
           arms: tuple[str, str] = ("baseline", "candidate")) -> list[tuple[dict, dict]]:
    """`count` slots of one pair each; in v2 a slot is retried until clean."""
    pairs = []
    for i in range(1, count + 1):
        for k in range(RETRIES + 1):
            stem = f"{tag}r{i}" + (f"-retry{k}" if k else "")
            pa, pb = d / f"{stem}-{arms[0]}.price.json", d / f"{stem}-{arms[1]}.price.json"
            a = price(arms[0], p, pa, extra, env)
            b = price(arms[1], p, pb, extra, env)
            if not pinned(a, b):
                (d / "PIN-MISMATCH").write_text(f"{stem}\n")
                raise SystemExit(f"pin mismatch at {d} {stem}: stop, the candidate is wrong")
            if clean(pa) and clean(pb):
                pairs.append((a, b))
                break
            print(f"{d.relative_to(RUNS)} {stem}: contended, slot run again", flush=True)
    return pairs


def main_comparison() -> None:
    for a, n in SIZES:
        d0 = RUNS / "main" / f"k{a}n{n}"
        d0.mkdir(parents=True, exist_ok=True)
        (d0 / "uptime-before.txt").exists() or (d0 / "uptime-before.txt").write_text(sh(["uptime"]) + "\n")
        for j in SETS:
            d = d0 / f"M{j}"
            p = params(a, n, j)
            pairs = rounds(d, p, "", [], ROUNDS)
            totals = [x["median"]["total_units"] for x, _ in pairs]
            spread = max(totals) / min(totals) if totals else float("inf")
            if spread > SPREAD_LIMIT:
                rounds(d, p, "double-", ["--repeats", "6", "--repeats-fast", "30"], ROUNDS)
            print(f"k{a}n{n} M{j}: A/A spread {spread:.3f}", flush=True)
        (d0 / "uptime-after.txt").exists() or (d0 / "uptime-after.txt").write_text(sh(["uptime"]) + "\n")


def aa() -> None:
    """v2: the noise floor, the baseline against a byte-identical copy."""
    if sha256(binary(ARMS["copy"])) != sha256(binary(ARMS["baseline"])):
        raise SystemExit("IC_BASELINE_COPY is not byte-identical to IC_BASELINE")
    for a, n in SIZES:
        d = RUNS / "aa" / f"k{a}n{n}" / "M1"
        rounds(d, params(a, n, 1), "", [], ROUNDS, ONE, ("baseline", "copy"))
        print(f"k{a}n{n} M1: A/A done", flush=True)


def control1() -> None:
    d = RUNS / "control1"
    d.mkdir(parents=True, exist_ok=True)
    for a, n in SIZES:
        p = params(a, n, 1)
        label = f"k{a}n{n}-M1"
        verdict = d / f"{label}-control.json"
        if verdict.exists():
            continue
        report = RUNS / "main" / f"k{a}n{n}" / "M1" / "r1-candidate.price.json"
        run_dir = d / f"{label}-wf"
        shutil.rmtree(run_dir, ignore_errors=True)
        wf = d / f"{label}-workflow.json"
        environ, cpus = ONE
        subprocess.run(["taskset", "-c", cpus, str(binary(ARMS["candidate"])), "workflow", "--params", str(p),
                        "--dir", str(run_dir), "--json", "--out", str(wf)], env=environ,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        res = subprocess.run([sys.executable, str(S20 / "control.py"), str(report), str(wf), str(run_dir)],
                             capture_output=True, text=True, check=False)
        verdict.write_text(res.stdout)
        shutil.rmtree(run_dir, ignore_errors=True)
        print(label, json.loads(res.stdout)["pass"], flush=True)


def rho() -> None:
    """§20 priced M1 with these seeds; the baseline must count the same."""
    d = RUNS / "rho"
    seeds = ["--rho-seed", str(0x200000 + 201), "--cold-rho-seed", str(0x210000)]
    out = d / "k0n41-M1-baseline.price.json"
    for k in range(RETRIES + 1):
        if k:
            out = d / f"k0n41-M1-retry{k}-baseline.price.json"
        rep = price("baseline", params(0, 41, 1), out, seeds)
        if clean(out):
            break
    s20 = json.loads((S20 / "runs" / "k0n41" / "measure" / "M1.price.json").read_text())
    keys = ("gae", "setup_gae", "steps_per_target", "counters", "solved_by", "recovered", "all_verified",
            "agrees_with_expected")
    same = {k: rep["rho_batch"].get(k) == s20["rho_batch"].get(k) for k in keys}
    doc = {"control": "baseline batch rho counts == §20's at k0n41 M1", "pass": all(same.values()),
           "report": out.name, "clean": clean(out), "compared": same,
           "s_priced_per_target": {"s20": s20["rho_batch"].get("s_priced_per_target"),
                                   "baseline": rep["rho_batch"].get("s_priced_per_target")},
           "pipeline_counts_equal_s20": rep.get("counts") == s20.get("counts"),
           "pipeline_recovered_equal_s20": rep.get("recovered") == s20.get("recovered")}
    (d / "k0n41-M1-verdict.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))


def threads() -> None:
    if ISOLATE:
        d, env, count = RUNS / "threads" / "k0n41-M1-t3", THREE, "three"
    else:
        d, env, count = RUNS / "threads" / "k0n41-M1", FOUR, "four"
    rounds(d, params(0, 41, 1), "", ["--allow-threads"], 3, env)
    print(f"k0n41 M1: {count} threads done", flush=True)


def probe() -> None:
    """The descent probe on both libraries, baseline then candidate, with
    the log tables committed at the declaration."""
    d = RUNS / "probe"
    d.mkdir(parents=True, exist_ok=True)
    for a, n in PROBED:
        logs = HERE / "probe" / f"k{a}n{n}-M1.logs.json"
        for arm in ("baseline", "candidate"):
            for k in range(RETRIES + 1):
                out = d / (f"k{a}n{n}-M1-{arm}.json" if not k else f"k{a}n{n}-M1-retry{k}-{arm}.json")
                if not out.exists():
                    cmd = [str(binary(PROBES[arm])), str(params(a, n, 1)), str(logs), "folded", "5"]
                    with open(out, "w") as stdout:
                        launch(cmd, ONE, out, stdout, out.with_suffix(".stderr"))
                    print(out.name, flush=True)
                if clean(out):
                    break


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    steps = {"manifest": manifest, "inputs": inputs, "aa": aa, "main": main_comparison, "control1": control1,
             "rho": rho, "threads": threads, "probe": probe}
    order = (("manifest", "inputs", "aa", "main", "rho", "threads", "probe") if ISOLATE
             else ("manifest", "inputs", "main", "control1", "rho", "threads", "probe"))
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "all":
        for name in order:
            print(f"== {name}", flush=True)
            steps[name]()
    elif cmd in steps:
        steps[cmd]()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
