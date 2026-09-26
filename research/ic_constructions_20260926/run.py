#!/usr/bin/env python3
"""Ledger §21's runs, in the order PROTOCOL.md declares them.

    python3 run.py manifest        # host manifest -> host.json
    python3 run.py inputs          # check §20's parameter files against inputs.sha256
    python3 run.py main            # 36 files x 5 rounds, baseline then candidate
    python3 run.py control1        # ic workflow == ic price, candidate, M1 at every size
    python3 run.py rho             # candidate re-prices batch rho at n = 41, M1
    python3 run.py threads         # 4 threads, n = 53 and n = 61, M1, 3 ABAB rounds
    python3 run.py constructions   # per-constructor prices, both binaries, n = 41 and 61
    python3 run.py all             # every step, in that order

The two binaries are named by the environment: IC_BASELINE and
IC_CANDIDATE for `ic`, PRICES_BASELINE and PRICES_CANDIDATE for
`examples/koblitz_construction_prices`.  Neither is in the tree; the
manifest records each one's sha256 and the commit it was built from.

Every single-thread process runs with RAYON_NUM_THREADS=1 under
`taskset -c 2`, one at a time.  Outputs go to runs/, one file per
process, and an existing file is never overwritten: a rerun resumes
where the last one stopped.  Nothing is dropped; a failed run keeps its
report and stderr.  A pin mismatch (counts or recovered logarithms
differing between the arms of a pair) stops the run, as declared.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = HERE / "runs"
S20 = ROOT / "research" / "ic_exponent_20260926"
# §20's nine sizes, in its order (by the size of r).
SIZES = [(1, 19), (1, 23), (1, 45), (0, 37), (1, 43), (1, 47), (0, 41), (0, 53), (0, 61)]
SETS = (1, 2, 3, 4)
ROUNDS = 5
SPREAD_LIMIT = 1.25
ARMS = {"baseline": "IC_BASELINE", "candidate": "IC_CANDIDATE"}
PRICES = {"baseline": "PRICES_BASELINE", "candidate": "PRICES_CANDIDATE"}
ONE = ({**os.environ, "RAYON_NUM_THREADS": "1"}, ["taskset", "-c", "2"])
FOUR = ({**os.environ, "RAYON_NUM_THREADS": "4"}, ["taskset", "-c", "0-3"])


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


def manifest() -> None:
    out = HERE / "host.json"
    if out.exists():
        print(f"{out} exists; not overwritten")
        return
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
                  "prices_sha256": sha256(binary(PRICES[arm])),
                  "built_from": os.environ.get(f"{ARMS[arm]}_COMMIT")}
            for arm in ARMS
        },
        "build": "cargo build --release --bin ic --example koblitz_construction_prices, clean tree, rustc above",
        "uptime": sh(["uptime"]),
        "pinning": "RAYON_NUM_THREADS=1 under taskset -c 2; the thread check RAYON_NUM_THREADS=4 under taskset -c 0-3",
        "hardware_class": "one x86-64 cloud container; no claim for Arm64, GPUs or other hosts",
    }
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
        environ, pin = env
        with open(out.with_suffix(".stderr"), "w") as e:
            subprocess.run([*pin, str(binary(ARMS[arm])), "price", "--params", str(p), "--json",
                            "--out", str(out), *extra], env=environ, stdout=subprocess.DEVNULL, stderr=e,
                           check=False)
    return json.loads(out.read_text())


def pinned(a: dict, b: dict) -> bool:
    return (a.get("status") == b.get("status") == "complete"
            and a.get("counts") == b.get("counts") and a.get("recovered") == b.get("recovered")
            and a.get("all_verified") is True and b.get("all_verified") is True)


def rounds(d: Path, p: Path, tag: str, extra: list[str], count: int, env=ONE) -> list[tuple[dict, dict]]:
    pairs = []
    for i in range(1, count + 1):
        a = price("baseline", p, d / f"{tag}r{i}-baseline.price.json", extra, env)
        b = price("candidate", p, d / f"{tag}r{i}-candidate.price.json", extra, env)
        if not pinned(a, b):
            (d / "PIN-MISMATCH").write_text(f"{tag}r{i}\n")
            raise SystemExit(f"pin mismatch at {d} {tag}r{i}: stop, the candidate is wrong")
        pairs.append((a, b))
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
            spread = max(totals) / min(totals)
            if spread > SPREAD_LIMIT:
                rounds(d, p, "double-", ["--repeats", "6", "--repeats-fast", "30"], ROUNDS)
            print(f"k{a}n{n} M{j}: A/A spread {spread:.3f}", flush=True)
        (d0 / "uptime-after.txt").exists() or (d0 / "uptime-after.txt").write_text(sh(["uptime"]) + "\n")


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
        environ, pin = ONE
        subprocess.run([*pin, str(binary(ARMS["candidate"])), "workflow", "--params", str(p), "--dir",
                        str(run_dir), "--json", "--out", str(wf)], env=environ, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False)
        res = subprocess.run([sys.executable, str(S20 / "control.py"), str(report), str(wf), str(run_dir)],
                             capture_output=True, text=True, check=False)
        verdict.write_text(res.stdout)
        shutil.rmtree(run_dir, ignore_errors=True)
        print(label, json.loads(res.stdout)["pass"], flush=True)


def rho() -> None:
    """§20 priced M1 with these seeds; the candidate must count the same."""
    d = RUNS / "rho"
    out = d / "k0n41-M1-candidate.price.json"
    rep = price("candidate", params(0, 41, 1), out, ["--rho-seed", str(0x200000 + 201),
                                                     "--cold-rho-seed", str(0x210000)])
    s20 = json.loads((S20 / "runs" / "k0n41" / "measure" / "M1.price.json").read_text())
    keys = ("gae", "setup_gae", "steps_per_target", "counters", "solved_by", "recovered", "all_verified",
            "agrees_with_expected")
    same = {k: rep["rho_batch"].get(k) == s20["rho_batch"].get(k) for k in keys}
    doc = {"control": "candidate batch rho counts == §20's at k0n41 M1", "pass": all(same.values()),
           "compared": same,
           "s_priced_per_target": {"s20": s20["rho_batch"].get("s_priced_per_target"),
                                   "candidate": rep["rho_batch"].get("s_priced_per_target")},
           "pipeline_counts_equal_s20": rep.get("counts") == s20.get("counts"),
           "pipeline_recovered_equal_s20": rep.get("recovered") == s20.get("recovered")}
    (d / "k0n41-M1-verdict.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))


def threads() -> None:
    for a, n in [(0, 53), (0, 61)]:
        d = RUNS / "threads" / f"k{a}n{n}-M1"
        rounds(d, params(a, n, 1), "", ["--allow-threads"], 3, FOUR)
        print(f"k{a}n{n} M1: four threads done", flush=True)


def constructions() -> None:
    d = RUNS / "constructions"
    d.mkdir(parents=True, exist_ok=True)
    environ, pin = ONE
    for a, n in [(0, 41), (0, 61)]:
        for arm in ("baseline", "candidate"):
            out = d / f"k{a}n{n}-M1-{arm}.json"
            if out.exists():
                continue
            res = subprocess.run([*pin, str(binary(PRICES[arm])), str(params(a, n, 1)), "5"], env=environ,
                                 capture_output=True, text=True, check=False)
            out.write_text(res.stdout)
            out.with_suffix(".stderr").write_text(res.stderr)
            print(out.name, flush=True)


def main() -> None:
    RUNS.mkdir(exist_ok=True)
    steps = {"manifest": manifest, "inputs": inputs, "main": main_comparison, "control1": control1,
             "rho": rho, "threads": threads, "constructions": constructions}
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "all":
        for name in ("manifest", "inputs", "main", "control1", "rho", "threads", "constructions"):
            print(f"== {name}", flush=True)
            steps[name]()
    elif cmd in steps:
        steps[cmd]()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
