#!/usr/bin/env python3
"""Ledger §20's runs, in the order PROTOCOL.md declares them.

    python3 run.py manifest            # host manifest -> host.json
    python3 run.py sweep <a> <n>       # the column / descent sweep on seed set W
    python3 run.py measure <a> <n>     # M1-M4 at the swept recipe, rho, controls
    python3 run.py controls            # the frozen headline and the thread's recipes
    python3 run.py all                 # every size, in the declared order

Every `ic price` runs with RAYON_NUM_THREADS=1 under `taskset -c 2`;
every control `ic workflow` runs the same way, after it, never beside
it.  Outputs go to runs/, one file per run, and an existing file is never
overwritten: a rerun resumes where the last one stopped.  Nothing is
dropped; a failed run keeps its report and stderr.
"""
from __future__ import annotations

import json
import math
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import make_params

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = HERE / "runs"
IC = ROOT / "target" / "release" / "ic"
PREDICTION = json.loads((HERE / "prediction.json").read_text())
SIZES = [(row["a"], row["n"]) for row in PREDICTION["rows"]]
SWEEP_SEED = 101
MEASURE_SEEDS = [201, 202, 203, 204]
SPREAD_LIMIT = 1.25
ENV = {**os.environ, "RAYON_NUM_THREADS": "1"}
PIN = ["taskset", "-c", "2"]


def sh(cmd: list[str]) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout.strip()


def manifest() -> None:
    out = HERE / "host.json"
    cpu = {}
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if ":" in line:
            k, v = (s.strip() for s in line.split(":", 1))
            if k in ("model name", "flags") and k not in cpu:
                cpu[k] = v
    flags = cpu.get("flags", "").split()
    wanted = [f for f in flags if f in ("popcnt", "avx2", "pclmulqdq", "bmi2", "gfni", "vpclmulqdq") or f.startswith("avx512")]
    mem = next((l for l in Path("/proc/meminfo").read_text().splitlines() if l.startswith("MemTotal")), "")
    doc = {
        "commit": sh(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        "tree_clean": sh(["git", "-C", str(ROOT), "status", "--porcelain"]) == "",
        "rustc": sh(["rustc", "--version"]),
        "cpu_model": cpu.get("model name"),
        "cpu_flags_relevant": wanted,
        "logical_cores": os.cpu_count(),
        "memory": mem,
        "os": platform.platform(),
        "arch": platform.machine(),
        "ic_binary_blake3": sh(["b3sum", "--no-names", str(IC)]) or None,
        "ic_binary_sha256": sh(["sha256sum", str(IC)]).split(" ")[0],
        "uptime": sh(["uptime"]),
        "pinning": "RAYON_NUM_THREADS=1, taskset -c 2",
        "hardware_class": "one x86-64 cloud container; no claim for Arm64, GPUs or other hosts",
    }
    if out.exists():
        print(f"{out} exists; not overwritten")
        return
    out.write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))


def price(params: Path, out: Path, extra: list[str]) -> dict:
    """Run `ic price` once into `out` unless it exists; return the report."""
    if not out.exists():
        err = out.with_suffix(".stderr")
        with open(err, "w") as e:
            subprocess.run([*PIN, str(IC), "price", "--params", str(params), "--json", "--out", str(out), *extra],
                           env=ENV, stdout=subprocess.DEVNULL, stderr=e, check=False)
    return json.loads(out.read_text())


def price_with_spread_rule(params: Path, out: Path, extra: list[str]) -> dict:
    """The declared spread rule: above 1.25, rerun once with double the
    repetitions and keep both; the rerun is the figure."""
    rep = price(params, out, extra)
    if rep.get("status") == "complete" and rep.get("spread_max_over_min", 0) > SPREAD_LIMIT:
        again = out.with_name(out.stem + "-double.json")
        rep = price(params, again, [*extra, "--repeats", "6", "--repeats-fast", "30"])
    return rep


def workflow_control(params: Path, out_dir: Path, price_report: Path) -> dict:
    verdict = out_dir / (params.stem + "-control.json")
    if verdict.exists():
        return json.loads(verdict.read_text())
    run_dir = out_dir / (params.stem + "-wf")
    shutil.rmtree(run_dir, ignore_errors=True)
    wf_report = out_dir / (params.stem + "-workflow.json")
    if wf_report.exists():
        wf_report.unlink()
    subprocess.run([*PIN, str(IC), "workflow", "--params", str(params), "--dir", str(run_dir), "--json",
                    "--out", str(wf_report)], env=ENV, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    res = subprocess.run([sys.executable, str(HERE / "control.py"), str(price_report), str(wf_report), str(run_dir)],
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
    """The base the thread's builder makes for a request (Amendment 1)."""
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
            rep = price(params, d / f"c{c}-m{m}.price.json", [])
            results[(c, m)] = rep

    for c in grid:
        run_point(c)
    while True:
        ok = {k: v for k, v in results.items() if v.get("status") == "complete"}
        if not ok:
            break
        best = min(ok, key=lambda k: ok[k]["median"]["s_per_target"])
        cols = sorted({k[0] for k in results})
        if extensions >= 3:
            break
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
                    "points": v.get("counts", {}).get("select", {}).get("points") if isinstance(v.get("counts"), dict) else None}
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
        out = d / f"M{seed - 200}.price.json"
        price_with_spread_rule(params, out, extra)
        workflow_control(params, d, out)
    (d / "uptime-after.txt").exists() or (d / "uptime-after.txt").write_text(sh(["uptime"]) + "\n")


def controls() -> None:
    d = RUNS / "controls"
    d.mkdir(parents=True, exist_ok=True)
    # Control 2: the frozen headline as it stands (its old rho baseline
    # removed from the workflow's copy only, so the control does not walk
    # the implemented rho 32 times; the pricer never runs that baseline).
    frozen = json.loads((ROOT / "docs/ic/params/k0n41-least-on-u150.json").read_text())
    frozen.pop("baseline", None)
    p = d / "frozen-headline.params.json"
    if not p.exists():
        p.write_text(json.dumps(frozen, indent=1) + "\n")
    out = d / "frozen-headline.price.json"
    price_with_spread_rule(p, out, ["--rho-seed", str(0x200001)])
    workflow_control(p, d, out)
    # Control 3: the thread's own n = 41 and n = 53 recipes on M1-M4.
    recipes = {
        "thread-n41": "docs/ic/params/k0n41-least-on-u150.json",
        "thread-n53": "docs/ic/params/k0n53-subgroup-aimed.json",
    }
    for label, src in recipes.items():
        base = json.loads((ROOT / src).read_text())
        base.pop("baseline", None)
        for seed in MEASURE_SEEDS:
            q = dict(base)
            q["name"] = f"{label}-s{seed}"
            q["seed"] = seed
            q["factor_base"] = {"mode": "spec", "spec": {**base["factor_base"]["spec"], "seed": seed}}
            q["targets"] = [{"random_seed": 100 * seed + i} for i in range(32)]
            p = d / f"{label}-M{seed - 200}.params.json"
            if not p.exists():
                p.write_text(json.dumps(q, indent=1) + "\n")
            out = d / f"{label}-M{seed - 200}.price.json"
            price_with_spread_rule(p, out, ["--rho-seed", str(0x200000 + seed)])
            workflow_control(p, d, out)


def main() -> None:
    RUNS.mkdir(exist_ok=True)
    cmd = sys.argv[1]
    if cmd == "manifest":
        manifest()
    elif cmd == "sweep":
        print(json.dumps(sweep(int(sys.argv[2]), int(sys.argv[3]))["chosen"]))
    elif cmd == "measure":
        measure(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "controls":
        controls()
    elif cmd == "all":
        manifest()
        for a, n in SIZES:
            print(f"k{a}n{n}: sweep", flush=True)
            print(json.dumps(sweep(a, n)["chosen"]), flush=True)
            print(f"k{a}n{n}: measure", flush=True)
            measure(a, n)
        print("controls", flush=True)
        controls()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
