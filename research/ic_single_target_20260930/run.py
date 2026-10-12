#!/usr/bin/env python3
"""Ledger §23's runs, in the order PROTOCOL.md declares them.

    python3 run.py manifest        # host manifest -> <runs>/host.json
    python3 run.py pin             # the batch pricer on §20's M1 files against §22's record
    python3 run.py size <a> <n>    # T01-T64 once, then T01-T04 again (the A/A)
    python3 run.py all             # manifest, pin, then every size in order of r

Every timed process is one `ic price --single-target` on one target
file, run through `tools/isolated_bench.py run --wait --cpus 2` with
`RAYON_NUM_THREADS=1`, after PSI `some avg10` has fallen below 4.0.  A
refused start is logged and tried again after 15 s.  A contended or
failed run is kept and run again up to twice; the first clean run is the
row.  The pin is counts only and runs untimed under `taskset -c 2`.

The binary is IC (IC_COMMIT names its commit); outputs go to IC_RUNS
(default runs/).  An existing file is never overwritten: a rerun resumes
where the last one stopped, and nothing is dropped.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import make_params  # noqa: E402

S20 = ROOT / "research" / "ic_exponent_20260926" / "runs"
S22 = ROOT / "research" / "ic_descent_20260930" / "runs-isolated" / "main"
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
TOOL = ROOT / "tools" / "isolated_bench.py"
TARGETS = 64
AA_TARGETS = 4
RHO_SEED_BASE = 0x230000
RETRIES = 2
REFUSAL_WAIT_S, REFUSAL_LIMIT = 15, 60
PSI_READY, PSI_WAIT_MAX_S = 4.0, 120
CPUS = "2"
ENV = {**os.environ, "RAYON_NUM_THREADS": "1"}
PIN_SIZES = [(1, 47), (0, 41), (0, 53), (0, 61)]


def sh(cmd: list[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout.strip()
    except FileNotFoundError:
        return ""


def ic() -> Path:
    path = os.environ.get("IC")
    if not path or not Path(path).exists():
        raise SystemExit("set IC to the round's ic binary")
    return Path(path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def psi_some_avg10() -> float:
    worst = 0.0
    for kind in ("cpu", "memory"):
        try:
            for line in Path(f"/proc/pressure/{kind}").read_text().splitlines():
                if line.startswith("some"):
                    worst = max(worst, float(line.split("avg10=")[1].split()[0]))
        except OSError:
            pass
    return worst


def wait_for_quiet() -> float:
    """Wait until PSI some avg10 is below PSI_READY; return the seconds waited."""
    start = time.monotonic()
    while psi_some_avg10() >= PSI_READY and time.monotonic() - start < PSI_WAIT_MAX_S:
        time.sleep(1)
    return time.monotonic() - start


def record_path(out: Path) -> Path:
    return out.with_name(out.name[: -len(".price.json")] + ".isolation.jsonl")


def clean(out: Path) -> bool:
    rec = record_path(out)
    if not rec.exists():
        return False
    run = json.loads(rec.read_text().splitlines()[-1])["run"]
    return run["exit_status"] == 0 and not run["contended"]


def launch(cmd: list[str], out: Path) -> None:
    """One timed process through the isolation tool, after PSI has fallen."""
    rec, label = record_path(out), str(out.relative_to(RUNS))
    stderr = out.with_name(out.name[: -len(".price.json")] + ".stderr")
    for attempt in range(REFUSAL_LIMIT):
        waited = wait_for_quiet()
        tmp = stderr.with_suffix(".stderr.tmp")
        with open(tmp, "w") as e:
            subprocess.run([sys.executable, str(TOOL), "run", "--wait", "--cpus", CPUS, "--out", str(rec),
                            "--label", label, "--", *cmd], env=ENV, stdout=subprocess.DEVNULL, stderr=e, check=False)
        if rec.exists():
            tmp.replace(stderr)
            with open(RUNS / "psi-waits.log", "a") as log:
                log.write(f"{label} {waited:.1f}\n")
            return
        with open(RUNS / "refusals.log", "a") as log:
            log.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {label} attempt {attempt + 1} "
                      f"(after {waited:.1f} s of PSI wait): {tmp.read_text().strip()}\n")
        tmp.unlink()
        time.sleep(REFUSAL_WAIT_S)
    raise SystemExit(f"{label}: the machine stayed busy for {REFUSAL_LIMIT} tries; stopped (resumable)")


def price_once(params: Path, out: Path, extra: list[str]) -> dict:
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        launch([str(ic()), "price", "--params", str(params), "--json", "--out", str(out), *extra], out)
    try:
        return json.loads(out.read_text())
    except (OSError, json.JSONDecodeError):
        return {"status": "no report"}


def price(params: Path, out: Path, extra: list[str]) -> dict:
    """The first clean run of up to 1 + RETRIES; every attempt is kept."""
    attempt = out
    rep = price_once(params, attempt, extra)
    for k in range(1, RETRIES + 1):
        if clean(attempt) and rep.get("status") == "complete":
            return rep
        attempt = out.with_name(out.name[: -len(".price.json")] + f"-retry{k}.price.json")
        rep = price_once(params, attempt, extra)
    return rep


def size_dir(a: int, n: int) -> Path:
    return RUNS / f"k{a}n{n}"


def target_file(a: int, n: int, i: int) -> Path:
    path = size_dir(a, n) / f"T{i:02d}.params.json"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(make_params.params(a, n, i), indent=1) + "\n")
    return path


def one(a: int, n: int, i: int, run: int) -> dict:
    extra = ["--single-target", "--rho-seed", str(RHO_SEED_BASE + i)]
    return price(target_file(a, n, i), size_dir(a, n) / f"T{i:02d}-R{run}.price.json", extra)


def size(a: int, n: int) -> None:
    d = size_dir(a, n)
    d.mkdir(parents=True, exist_ok=True)
    stop = d / "STOPPED"
    if stop.exists():
        print(f"k{a}n{n}: stopped earlier ({stop.read_text().strip()}); not resumed")
        return
    (d / "uptime-before.txt").exists() or (d / "uptime-before.txt").write_text(sh(["uptime"]) + "\n")
    order = [(i, 1) for i in range(1, TARGETS + 1)] + [(i, 2) for i in range(1, AA_TARGETS + 1)]
    for i, run in order:
        rep = one(a, n, i, run)
        if rep.get("status") != "complete":
            # PROTOCOL.md: a verification failure stops the size.
            stop.write_text(f"T{i:02d}-R{run}: status {rep.get('status')}\n")
            print(f"k{a}n{n}: T{i:02d}-R{run} is {rep.get('status')}; size stopped", flush=True)
            return
        print(f"k{a}n{n}: T{i:02d}-R{run} online speedup {rep['median']['online_speedup']:.2f}", flush=True)
    (d / "uptime-after.txt").exists() or (d / "uptime-after.txt").write_text(sh(["uptime"]) + "\n")


def pin() -> dict:
    """The batch pricer on §20's frozen M1 files: counts and recovered
    logarithms must equal §22's isolated record, or the pricer is not the
    one §22 measured."""
    d = RUNS / "pin"
    d.mkdir(parents=True, exist_ok=True)
    verdict_path = d / "pin.json"
    if verdict_path.exists():
        return json.loads(verdict_path.read_text())
    rows = []
    for a, n in PIN_SIZES:
        out = d / f"k{a}n{n}-M1.price.json"
        if not out.exists():
            subprocess.run(["taskset", "-c", CPUS, str(ic()), "price", "--params",
                            str(S20 / f"k{a}n{n}" / "measure" / "M1.params.json"), "--repeats", "1",
                            "--repeats-fast", "1", "--json", "--out", str(out)],
                           env=ENV, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        new = json.loads(out.read_text())
        old = json.loads((S22 / f"k{a}n{n}" / "M1" / "r1-candidate.price.json").read_text())
        rows.append({"curve": f"k{a}n{n}", "counts_equal": new["counts"] == old["counts"],
                     "recovered_equal": new["recovered"] == old["recovered"], "status": new["status"]})
    verdict = {"rows": rows, "held": all(r["counts_equal"] and r["recovered_equal"] for r in rows)}
    verdict_path.write_text(json.dumps(verdict, indent=1) + "\n")
    return verdict


def manifest() -> None:
    out = RUNS / "host.json"
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
    mem = next((line for line in Path("/proc/meminfo").read_text().splitlines() if line.startswith("MemTotal")), "")
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
        "ic_binary_sha256": sha256(ic()),
        "ic_built_from": os.environ.get("IC_COMMIT"),
        "isolation_tool": {"path": "tools/isolated_bench.py", "sha256": sha256(TOOL)},
        "pinning": ("tools/isolated_bench.py run --wait --cpus 2 with RAYON_NUM_THREADS=1, after PSI some avg10 "
                    "falls below 4.0; the tool's defaults (settle 2 s, other processes at most 0.10 CPUs, PSI at "
                    "most 5); the pin untimed under taskset -c 2"),
        "smt_siblings": {c: Path(f"/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list").read_text().strip()
                         for c in range(os.cpu_count() or 0)},
        "uptime": sh(["uptime"]),
        "hardware_class": "one x86-64 cloud container; no claim for Arm64, GPUs or other hosts",
    }
    out.write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "manifest":
        manifest()
    elif cmd == "pin":
        print(json.dumps(pin(), indent=1))
    elif cmd == "size":
        size(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "all":
        manifest()
        verdict = pin()
        if not verdict["held"]:
            raise SystemExit("the pin failed: the pricer is not the one §22 measured; stopped")
        for a, n in make_params.SIZES:
            print(f"k{a}n{n}: start", flush=True)
            size(a, n)
            print(f"k{a}n{n}: done", flush=True)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
