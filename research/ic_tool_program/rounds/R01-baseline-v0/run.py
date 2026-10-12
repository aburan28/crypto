#!/usr/bin/env python3
"""R01's steps, in the order PROTOCOL.md declares them.

    IC=<v0 binary> IC_COMMIT=<build commit> python3 run.py all
    python3 run.py <step>     # manifest | check | pin | smoke | profile | aa | thp | calib | callgrind

Outputs go to runs/ (or IC_RUNS).  Nothing is overwritten: every step
resumes where it stopped.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402

ROOT = bench.ROOT
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
S23 = ROOT / "research" / "ic_single_target_20260930" / "runs"
S23_SIZES = ["k1n47", "k0n57", "k0n41", "k0n53", "k1n59", "k0n61"]
THP_SIZES = ["k0n53", "k1n59", "k0n61"]
CALLGRIND_ROWS = ["k0n41-M1-T01", "k0n61-M1-T01"]
ROUNDS = 5
CALIB = ROOT / "research" / "notes" / "index-calculus" / "cachegrind_n41_20260929_run" / "mem_calib.c"
THP_ENV = {**bench.ENV, "GLIBC_TUNABLES": "glibc.malloc.hugetlb=1"}


def ic() -> Path:
    path = os.environ.get("IC")
    if not path or not Path(path).exists():
        raise SystemExit("set IC to v0's ic binary")
    return Path(path)


def copy_of_ic() -> Path:
    """A byte-identical copy for the A/A, next to the original."""
    src = ic()
    dst = src.with_name(src.name + "-aa-copy")
    if not dst.exists():
        shutil.copy2(src, dst)
    if bench.sha256(dst) != bench.sha256(src):
        raise SystemExit("the A/A copy is not byte-identical")
    return dst


def manifest() -> None:
    doc = bench.host_manifest(RUNS / "host.json", {
        "v0": {"path_basename": ic().name, "sha256": bench.sha256(ic()),
               "built_from": os.environ.get("IC_COMMIT"), "src_tree": "003badc259bcef6e56ad19ac8ae93df4221ffe6e"},
    })
    print(json.dumps({k: doc[k] for k in ("cpu_model", "logical_cores", "transparent_hugepage", "binaries")},
                     indent=1))


def check() -> None:
    r = subprocess.run([sys.executable, str(bench.SUITE / "make_suite.py"), "--check"], capture_output=True,
                       text=True)
    if r.returncode != 0:
        raise SystemExit(f"suite v1 does not re-derive: {r.stdout}")


def rows_by_id() -> dict[str, dict]:
    return {r["id"]: r for r in json.loads((bench.SUITE / "SUITE.json").read_text())["rows"]}


def pin() -> dict:
    out = RUNS / "pin" / "pin.json"
    if out.exists():
        return json.loads(out.read_text())
    rows = rows_by_id()
    verdict = []
    for size in S23_SIZES:
        for t in (1, 2):
            row = rows[f"{size}-M1-T{t:02d}"]
            new = bench.untimed(ic(), row, RUNS / "pin" / f"{row['id']}.price.json")
            old = json.loads((S23 / size / f"T{t:02d}-R1.price.json").read_text())
            a, b = bench.outputs(new), bench.outputs(old)
            verdict.append({"row": row["id"], "equal": a == b,
                            "differs_in": sorted(k for k in a if a[k] != b[k])})
    doc = {"rows": verdict, "held": all(v["equal"] for v in verdict),
           "reference": "research/ic_single_target_20260930/runs/<size>/T0{1,2}-R1.price.json"}
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def smoke() -> dict:
    out = RUNS / "smoke" / "smoke.json"
    if out.exists():
        return json.loads(out.read_text())
    result = []
    for row in bench.suite_rows("smoke"):
        rep = bench.untimed(ic(), row, RUNS / "smoke" / f"{row['id']}.price.json")
        result.append({"row": row["id"], **bench.outputs(rep)})
    doc = {"rows": result, "ok": all(r["status"] == "complete" and r["all_verified"] for r in result)}
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def profile() -> None:
    bench.interleave({"v0": ic()}, bench.suite_rows("S"), 1, RUNS / "profile")


def m1_rows(sizes: list[str] | None = None) -> list[dict]:
    rows = [r for r in bench.suite_rows("S") if r["recipe_seed"] == 201]
    return [r for r in rows if sizes is None or f"k{r['a']}n{r['n']}" in sizes]


def aa() -> None:
    bench.interleave({"A": ic(), "A2": copy_of_ic()}, m1_rows(), ROUNDS, RUNS / "aa")


def thp() -> None:
    bench.interleave({"4k": ic(), "thp": ic()}, m1_rows(THP_SIZES), ROUNDS, RUNS / "thp",
                     envs={"4k": bench.ENV, "thp": THP_ENV})


def calib() -> None:
    d = RUNS / "calib"
    d.mkdir(parents=True, exist_ok=True)
    binary = d / "mem_calib"
    if not binary.exists():
        subprocess.run(["gcc", "-O2", "-o", str(binary), str(CALIB)], check=True)
    (d / "source.json").exists() or (d / "source.json").write_text(json.dumps(
        {"source": str(CALIB.relative_to(ROOT)), "sha256": bench.sha256(CALIB), "cc": bench.sh(["gcc", "--version"]).splitlines()[0],
         "flags": "-O2"}, indent=1) + "\n")
    kib, mib = 1024, 1024 * 1024
    sizes = [16 * kib, 256 * kib, mib, 2 * mib, 4 * mib, 8 * mib, 16 * mib, 32 * mib, 64 * mib, 256 * mib,
             1024 * mib]
    for mode, env in (("4k", bench.ENV), ("thp", THP_ENV)):
        for w in sizes:
            steps = 20_000_000 if w <= mib else 3_000_000 if w <= 64 * mib else 1_500_000
            out = d / mode / f"chase-{w}.json"
            if out.exists():
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            tmp = out.with_suffix(".stdout")
            bench.launch(["sh", "-c", f"{binary} chase {w} 1 {steps} 5 {steps // 10} > {tmp}"], out, RUNS, env)
            lines = tmp.read_text().strip().splitlines()
            out.write_text((lines[-1] if lines else "{}") + "\n")
        for k in (1, 4, 8, 16):
            out = d / mode / f"mlp-256MiB-k{k}.json"
            if out.exists():
                continue
            tmp = out.with_suffix(".stdout")
            bench.launch(["sh", "-c", f"{binary} chase {256 * mib} {k} 1000000 5 100000 > {tmp}"], out, RUNS, env)
            lines = tmp.read_text().strip().splitlines()
            out.write_text((lines[-1] if lines else "{}") + "\n")


def callgrind() -> None:
    d = RUNS / "callgrind"
    d.mkdir(parents=True, exist_ok=True)
    rows = rows_by_id()
    for rid in CALLGRIND_ROWS:
        row = rows[rid]
        cg = d / f"{rid}.callgrind.out"
        if cg.exists():
            continue
        rep = d / f"{rid}.price.json"
        cmd = ["valgrind", "--tool=callgrind", "--cache-sim=yes", "--I1=32768,8,64", "--D1=32768,8,64",
               "--LL=2097152,16,64", f"--callgrind-out-file={cg}",
               *bench.price_cmd(ic(), row, rep, ("--repeats", "1", "--repeats-fast", "1"))]
        with open(d / f"{rid}.valgrind.log", "w") as log:
            subprocess.run(["taskset", "-c", bench.CPUS, *cmd], env=bench.ENV, stdout=subprocess.DEVNULL,
                           stderr=log, check=False)
        with open(d / f"{rid}.annotate.txt", "w") as f:
            subprocess.run(["callgrind_annotate", "--inclusive=no", "--threshold=99", str(cg)], stdout=f,
                           stderr=subprocess.DEVNULL, check=False)
        with open(d / f"{rid}.annotate-inclusive.txt", "w") as f:
            subprocess.run(["callgrind_annotate", "--inclusive=yes", "--threshold=99", str(cg)], stdout=f,
                           stderr=subprocess.DEVNULL, check=False)


STEPS = {"manifest": manifest, "check": check, "pin": pin, "smoke": smoke, "profile": profile, "aa": aa,
         "thp": thp, "calib": calib, "callgrind": callgrind}


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd in STEPS:
        result = STEPS[cmd]()
        if isinstance(result, dict):
            print(json.dumps(result, indent=1))
        return
    if cmd != "all":
        raise SystemExit(__doc__)
    manifest()
    check()
    verdict = pin()
    print(json.dumps(verdict, indent=1), flush=True)
    if not verdict["held"]:
        raise SystemExit("the pin failed; R01 stops before the timed steps (PROTOCOL.md step 2)")
    if not smoke()["ok"]:
        raise SystemExit("a smoke row failed")
    for step in (profile, aa, thp, calib, callgrind):
        check()
        step()


if __name__ == "__main__":
    main()
