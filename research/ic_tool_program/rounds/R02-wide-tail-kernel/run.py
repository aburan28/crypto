#!/usr/bin/env python3
"""R02's steps, in the order PROTOCOL.md declares them.

    IC_BASE=<v0> IC_CAND=<candidate> IC_CAND_COMMIT=<commit> python3 run.py all
    python3 run.py <step>   # manifest | control | pin | holdouts | compare | holdout | callgrind | rule

Outputs go to runs/ (or IC_RUNS).  Nothing is overwritten: every step
resumes where it stopped.  `rule` runs only when the round is accepted.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402
import stats  # noqa: E402

ROOT = bench.ROOT
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
R01 = HERE.parent / "R01-baseline-v0" / "runs"
ROUNDS = 5
CONTROL_SIZE = "k0n53"
CALLGRIND_ROWS = ["k0n61-M1-T01", "k0n41-M1-T01"]
HOLDOUT_SEED = 205
HOLDOUT_TARGETS = (101, 102)
SCALAR_ENV = {**bench.ENV, "KIC_SCAN_SIMD": "0"}


def arm(name: str) -> Path:
    path = os.environ.get(name)
    if not path or not Path(path).exists():
        raise SystemExit(f"set {name}")
    return Path(path)


def manifest() -> None:
    doc = bench.host_manifest(RUNS / "host.json", {
        "v0": {"path_basename": arm("IC_BASE").name, "sha256": bench.sha256(arm("IC_BASE")),
               "built_from": "4afd29903e4fcf65c1d7096083bbe9b9f5ec0a66"},
        "candidate": {"path_basename": arm("IC_CAND").name, "sha256": bench.sha256(arm("IC_CAND")),
                      "built_from": os.environ.get("IC_CAND_COMMIT")},
    })
    r01 = json.loads((R01 / "host.json").read_text())
    same = all(doc[k] == r01[k] for k in ("cpu_model", "cpu_flags_relevant", "logical_cores", "memory",
                                          "transparent_hugepage", "os"))
    (RUNS / "aa-source.json").write_text(json.dumps(
        {"host_matches_r01": same, "aa": "R01's" if same else "R02's own (run `aa`)"}, indent=1) + "\n")
    print(json.dumps({"host_matches_r01": same}, indent=1))


def rows() -> list[dict]:
    return bench.suite_rows("S")


def m1(size: str) -> list[dict]:
    return [r for r in rows() if r["recipe_seed"] == 201 and f"k{r['a']}n{r['n']}" == size]


def control() -> dict:
    """KIC_SCAN_SIMD=0 against the default on v0 (PROTOCOL.md, "A control first")."""
    d = RUNS / "control"
    bench.interleave({"simd": arm("IC_BASE"), "scalar": arm("IC_BASE")}, m1(CONTROL_SIZE), ROUNDS, d,
                     envs={"simd": bench.ENV, "scalar": SCALAR_ENV})
    ratios = []
    for r in m1(CONTROL_SIZE):
        for k in range(1, ROUNDS + 1):
            a = bench.load(bench.figure_path(d / "scalar" / r["id"] / f"r{k}.price.json"))
            b = bench.load(bench.figure_path(d / "simd" / r["id"] / f"r{k}.price.json"))
            if a.get("status") == b.get("status") == "complete":
                ratios.append(stats.setup_phase_ns(a)["collect"] / stats.setup_phase_ns(b)["collect"])
    verdict = {"collect_scalar_over_simd": stats.geo_ci(ratios), "threshold": 1.3}
    verdict["confirmed"] = verdict["collect_scalar_over_simd"].get("geomean", 0) >= 1.3
    (d / "verdict.json").write_text(json.dumps(verdict, indent=1) + "\n")
    return verdict


def pin() -> dict:
    """The candidate's outputs on every suite row against v0's from R01."""
    out = RUNS / "pin" / "pin.json"
    if out.exists():
        return json.loads(out.read_text())
    result = []
    for row in bench.suite_rows("S") + bench.suite_rows("smoke"):
        new = bench.untimed(arm("IC_CAND"), row, RUNS / "pin" / f"{row['id']}.price.json")
        if row["tier"] == "S":
            old = bench.load(bench.figure_path(R01 / "profile" / "v0" / row["id"] / "r1.price.json"))
        else:
            old = bench.load(R01 / "smoke" / f"{row['id']}.price.json")
        a, b = bench.outputs(new), bench.outputs(old)
        result.append({"row": row["id"], "equal": a == b, "differs_in": sorted(k for k in a if a[k] != b[k])})
    doc = {"rows": result, "held": all(r["equal"] for r in result)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def holdout_rows() -> list[dict]:
    """Seed 205 and targets T101/T102 per size, by suite v1's own construction."""
    spec = importlib.util.spec_from_file_location("suite_v1", bench.SUITE / "make_suite.py")
    suite = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(suite)
    d = HERE / "holdouts"
    out = []
    for (a, n), (columns, m, _) in suite.RECIPES.items():
        for i in HOLDOUT_TARGETS:
            text = suite.text(suite.params(a, n, columns, m, HOLDOUT_SEED, i))
            path = d / f"k{a}n{n}" / f"M5-T{i}.json"
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
            elif path.read_text() != text:
                raise SystemExit(f"{path} differs from its construction")
            r = suite.s20.R[(a, n)]
            out.append({"id": f"k{a}n{n}-M5-T{i}", "tier": "holdout", "a": a, "n": n, "r": r,
                        "recipe_seed": HOLDOUT_SEED, "target": i, "rho_seed": suite.RHO_SEED_BASE + i,
                        # Absolute, so bench.price_cmd's join with the suite leaves it as it is.
                        "params": str(path.resolve())})
    return out


def compare() -> None:
    bench.interleave({"v0": arm("IC_BASE"), "cand": arm("IC_CAND")}, rows(), ROUNDS, RUNS / "compare")


def holdout() -> None:
    bench.interleave({"v0": arm("IC_BASE"), "cand": arm("IC_CAND")}, holdout_rows(), ROUNDS, RUNS / "holdout")


def callgrind() -> None:
    d = RUNS / "callgrind"
    d.mkdir(parents=True, exist_ok=True)
    by_id = {r["id"]: r for r in rows()}
    for name, binary in (("v0", arm("IC_BASE")), ("cand", arm("IC_CAND"))):
        for rid in CALLGRIND_ROWS:
            cg = d / f"{name}-{rid}.callgrind.out"
            if cg.exists():
                continue
            rep = d / f"{name}-{rid}.price.json"
            cmd = ["valgrind", "--tool=callgrind", "--cache-sim=yes", "--I1=32768,8,64", "--D1=32768,8,64",
                   "--LL=2097152,16,64", f"--callgrind-out-file={cg}",
                   *bench.price_cmd(binary, by_id[rid], rep, ("--repeats", "1", "--repeats-fast", "1"))]
            with open(d / f"{name}-{rid}.valgrind.log", "w") as log:
                subprocess.run(["taskset", "-c", bench.CPUS, *cmd], env=bench.ENV, stdout=subprocess.DEVNULL,
                               stderr=log, check=False)
            with open(d / f"{name}-{rid}.annotate.txt", "w") as f:
                subprocess.run(["callgrind_annotate", "--inclusive=no", "--threshold=99", str(cg)], stdout=f,
                               stderr=subprocess.DEVNULL, check=False)


def rule() -> None:
    """§23's protocol at its six sizes with 64 targets, on the candidate."""
    s23 = ROOT / "research" / "ic_single_target_20260930"
    env = {**os.environ, "IC": str(arm("IC_CAND")), "IC_COMMIT": os.environ.get("IC_CAND_COMMIT", ""),
           "IC_RUNS": str(RUNS / "rule")}
    for a, n in [(1, 47), (0, 57), (0, 41), (0, 53), (1, 59), (0, 61)]:
        subprocess.run([sys.executable, str(s23 / "run.py"), "size", str(a), str(n)], env=env, check=True)


STEPS = {"manifest": manifest, "control": control, "pin": pin, "holdouts": holdout_rows, "compare": compare,
         "holdout": holdout, "callgrind": callgrind, "rule": rule}


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd in STEPS:
        result = STEPS[cmd]()
        if isinstance(result, (dict, list)):
            print(json.dumps(result, indent=1))
        return
    if cmd != "all":
        raise SystemExit(__doc__)
    manifest()
    verdict = control()
    print(json.dumps(verdict, indent=1), flush=True)
    if not verdict["confirmed"]:
        raise SystemExit("the control did not confirm the mechanism; R02 stops (PROTOCOL.md)")
    if not pin()["held"]:
        raise SystemExit("an output differs; R02 stops (PROTOCOL.md)")
    compare()
    holdout()
    callgrind()


if __name__ == "__main__":
    main()
