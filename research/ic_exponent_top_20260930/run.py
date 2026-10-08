#!/usr/bin/env python3
"""Ledger §23's runs, in the order PROTOCOL.md declares them.

    python3 run.py manifest            # host manifest -> <runs>/host.json
    python3 run.py sweep <a> <n>       # the column / descent sweep on seed set W
    python3 run.py measure <a> <n>     # M1-M8 at the swept recipe, rho, Control 1 on M1
    python3 run.py all                 # every size, in the declared order

§20's rules and procedure (research/ic_exponent_20260926/make_params.py
and run.py), with two changes the protocol declares:

- every timed process (`ic price`, sweep and sets) runs through
  `tools/isolated_bench.py run --wait --cpus 2`, after the harness waits
  for PSI `some avg10` to fall below 4.0; a contended or failed run is
  kept and run again up to twice, and the first clean run is the figure;
- eight measurement sets, M1-M8 (seeds 201-208), and Control 1 on M1.

The binary is IC (with IC_COMMIT naming its commit); outputs go to
IC_RUNS (default runs/).  An existing file is never overwritten: a
rerun resumes where the last one stopped, and nothing is dropped.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S20 = ROOT / "research" / "ic_exponent_20260926"
sys.path.insert(0, str(S20))
import make_params  # noqa: E402  §20's rules, unchanged

# The two sizes §20 did not have, from curve_records.json.
make_params.R.update({(0, 57): 275295876199, (1, 59): 25179555920633})

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
TOOL = ROOT / "tools" / "isolated_bench.py"
PREDICTION = json.loads((HERE / "prediction.json").read_text())
SIZES = [(row["a"], row["n"]) for row in PREDICTION["rows"]]
SWEEP_SEED = 101
MEASURE_SEEDS = list(range(201, 209))
SPREAD_LIMIT = 1.25
RETRIES = 2
REFUSAL_WAIT_S, REFUSAL_LIMIT = 15, 60
PSI_READY, PSI_WAIT_MAX_S = 4.0, 120
CPUS = "2"
ENV = {**os.environ, "RAYON_NUM_THREADS": "1"}


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
    name = out.name
    for suffix in (".price.json", ".json"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return out.with_name(name + ".isolation.jsonl")


def clean(out: Path) -> bool:
    rec = record_path(out)
    if not rec.exists():
        return False
    run = json.loads(rec.read_text().splitlines()[-1])["run"]
    return run["exit_status"] == 0 and not run["contended"]


def launch(cmd: list[str], out: Path) -> None:
    """One timed process through the isolation tool, after PSI has fallen."""
    rec, label = record_path(out), str(out.relative_to(RUNS))
    stderr = out.with_suffix(".stderr")
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
    return json.loads(out.read_text())


def price(params: Path, out: Path, extra: list[str]) -> dict:
    """The first clean run of up to 1 + RETRIES; every attempt is kept."""
    base, attempt = out, out
    rep = price_once(params, attempt, extra)
    for k in range(1, RETRIES + 1):
        if clean(attempt):
            return rep
        attempt = base.with_name(base.name[: -len(".json")] + f"-retry{k}.json")
        rep = price_once(params, attempt, extra)
    return rep


def price_with_spread_rule(params: Path, out: Path, extra: list[str]) -> dict:
    """§20's rule: above 1.25, rerun once with double the repetitions and
    keep both; the rerun is the figure."""
    rep = price(params, out, extra)
    if rep.get("status") == "complete" and rep.get("spread_max_over_min", 0) > SPREAD_LIMIT:
        again = out.with_name(out.stem + "-double.json")
        rep = price(params, again, [*extra, "--repeats", "6", "--repeats-fast", "30"])
    return rep


def workflow_control(params: Path, out_dir: Path, price_report: Path) -> dict:
    """Control 1, counts only and untimed: after the pricing, never beside it."""
    verdict = out_dir / (params.stem + "-control.json")
    if verdict.exists():
        return json.loads(verdict.read_text())
    run_dir = out_dir / (params.stem + "-wf")
    shutil.rmtree(run_dir, ignore_errors=True)
    wf_report = out_dir / (params.stem + "-workflow.json")
    if wf_report.exists():
        wf_report.unlink()
    subprocess.run(["taskset", "-c", CPUS, str(ic()), "workflow", "--params", str(params), "--dir", str(run_dir),
                    "--json", "--out", str(wf_report)], env=ENV, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=False)
    res = subprocess.run([sys.executable, str(S20 / "control.py"), str(price_report), str(wf_report), str(run_dir)],
                         capture_output=True, text=True, check=False)
    verdict.write_text(res.stdout)
    shutil.rmtree(run_dir, ignore_errors=True)
    return json.loads(res.stdout)


def write_params(a: int, n: int, columns: int, m: int, seed: int, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(make_params.params(a, n, columns, m, seed), indent=1) + "\n")
    return path


def actual(columns: int) -> int:
    """The base the thread's builder makes for a request (§20, Amendment 1)."""
    return max(8, 8 * math.ceil(columns / 8))


def size_dir(a: int, n: int) -> Path:
    return RUNS / f"k{a}n{n}"


def sweep(a: int, n: int) -> dict:
    d = size_dir(a, n) / "sweep"
    d.mkdir(parents=True, exist_ok=True)
    chosen_path = d / "chosen.json"
    if chosen_path.exists():
        return json.loads(chosen_path.read_text())
    row = next(r for r in PREDICTION["rows"] if (r["a"], r["n"]) == (a, n))
    grid = sorted({actual(c) for c in row["sweep_grid_columns"]})
    results: dict[tuple[int, int], dict] = {}
    extensions = 0

    def run_point(c: int) -> None:
        for m in (2, 3):
            if (c, m) in results:
                continue
            params = write_params(a, n, c, m, SWEEP_SEED, d / f"c{c}-m{m}.params.json")
            results[(c, m)] = price(params, d / f"c{c}-m{m}.price.json", [])

    for c in grid:
        run_point(c)
    while True:
        ok = {k: v for k, v in results.items() if v.get("status") == "complete"}
        if not ok or extensions >= 3:
            break
        best = min(ok, key=lambda k: ok[k]["median"]["s_per_target"])
        cols = sorted({k[0] for k in results})
        if best[0] == cols[0] and cols[0] > 8:
            run_point(cols[0] - 8)
        elif best[0] == cols[-1]:
            run_point(actual(math.ceil(cols[-1] * math.sqrt(2))))
        else:
            break
        extensions += 1
    ok = {k: v for k, v in results.items() if v.get("status") == "complete"}
    best = min(ok, key=lambda k: ok[k]["median"]["s_per_target"])
    chosen = {
        "curve": f"k{a}n{n}", "seed_set": "W", "seed": SWEEP_SEED,
        "grid_requested": row["sweep_grid_columns"], "grid_actual": grid, "extensions": extensions,
        "points": [{"columns": c, "descent_summands": m, "status": v.get("status"), "message": v.get("message"),
                    "s_per_target": v.get("median", {}).get("s_per_target"),
                    "points": v.get("counts", {}).get("select", {}).get("points")
                    if isinstance(v.get("counts"), dict) else None}
                   for (c, m), v in sorted(results.items())],
        "chosen": {"columns": best[0], "descent_summands": best[1], "s_per_target": ok[best]["median"]["s_per_target"]},
    }
    chosen_path.write_text(json.dumps(chosen, indent=1) + "\n")
    return chosen


def measure(a: int, n: int) -> None:
    chosen = sweep(a, n)["chosen"]
    c, m = chosen["columns"], chosen["descent_summands"]
    d = size_dir(a, n) / "measure"
    d.mkdir(parents=True, exist_ok=True)
    (d / "uptime-before.txt").exists() or (d / "uptime-before.txt").write_text(sh(["uptime"]) + "\n")
    for seed in MEASURE_SEEDS:
        params = write_params(a, n, c, m, seed, d / f"M{seed - 200}.params.json")
        extra = ["--rho-seed", str(0x200000 + seed)]
        if seed == MEASURE_SEEDS[0]:
            extra += ["--cold-rho-seed", str(0x210000)]
        price_with_spread_rule(params, d / f"M{seed - 200}.price.json", extra)
    m1 = d / "M1.price.json"
    workflow_control(d / "M1.params.json", d, m1)
    (d / "uptime-after.txt").exists() or (d / "uptime-after.txt").write_text(sh(["uptime"]) + "\n")


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
        "ic_binary_sha256": sha256(ic()),
        "ic_built_from": os.environ.get("IC_COMMIT"),
        "isolation_tool": {"path": "tools/isolated_bench.py", "sha256": sha256(TOOL)},
        "pinning": ("tools/isolated_bench.py run --wait --cpus 2 with RAYON_NUM_THREADS=1, after PSI some avg10 "
                    "falls below 4.0; the tool's defaults (settle 2 s, other processes at most 0.10 CPUs, PSI at "
                    "most 5); Control 1 untimed under taskset -c 2"),
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
    elif cmd == "sweep":
        print(json.dumps(sweep(int(sys.argv[2]), int(sys.argv[3]))["chosen"]))
    elif cmd == "measure":
        measure(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "all":
        manifest()
        for a, n in SIZES:
            print(f"k{a}n{n}: sweep", flush=True)
            print(json.dumps(sweep(a, n)["chosen"]), flush=True)
            print(f"k{a}n{n}: measure", flush=True)
            measure(a, n)
            print(f"k{a}n{n}: done", flush=True)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
