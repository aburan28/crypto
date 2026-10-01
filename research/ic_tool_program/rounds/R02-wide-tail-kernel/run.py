#!/usr/bin/env python3
"""R02's steps, in the order PROTOCOL.md declares them (amendment 3).

    IC_V0=<R01's v0> IC_BASE=<v0'> IC_BASE_COMMIT=<commit> \\
        IC_CAND=<candidate> IC_CAND_COMMIT=<commit> python3 run.py all
    python3 run.py <step>

The steps are manifest, control, pin, compare, extend, holdout, v0check,
callgrind, aa and rule.  `aa` runs only when the host differs from R01's
(it needs IC_BASE_COPY, a byte-identical copy of IC_BASE), and `rule` only
when the round is accepted.

Outputs go to runs/ (or IC_RUNS), one directory per row, named
`<slug>/<recipe>-<target>` as AGENTS.md §11 names run files.  Nothing is
overwritten: every step resumes where it stopped.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402
import stats  # noqa: E402

ROOT = bench.ROOT
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
R01_DIR = HERE.parent / "R01-baseline-v0"
V0_COMMIT = "4afd29903e4fcf65c1d7096083bbe9b9f5ec0a66"
ROUNDS, EXTENDED_ROUNDS, HALF_WIDTH_LIMIT = 5, 10, 0.05
CONTROL = (0, 53)
CALLGRIND = [(0, 61), (1, 59), (0, 41)]
HOLDOUT_SEED = 205
HOLDOUT_TARGETS = (101, 102)
SCALAR_ENV = {**bench.ENV, "KIC_SCAN_SIMD": "0"}


def arm(name: str) -> Path:
    path = os.environ.get(name)
    if not path or not Path(path).exists():
        raise SystemExit(f"set {name}")
    return Path(path)


def r01_runs() -> Path:
    """R01's run tree, checked against its SHA-256 and extracted on first use."""
    runs = R01_DIR / "runs"
    if not runs.exists():
        want = (R01_DIR / "runs.tar.xz.sha256").read_text().split()[0]
        if bench.sha256(R01_DIR / "runs.tar.xz") != want:
            raise SystemExit("R01's runs.tar.xz does not match runs.tar.xz.sha256")
        with tarfile.open(R01_DIR / "runs.tar.xz") as tar:
            tar.extractall(R01_DIR, filter="data")
    return runs


def rows() -> list[dict]:
    return bench.slug_rows(bench.suite_rows("S"))


def m1(a: int | None = None, n: int | None = None) -> list[dict]:
    return [r for r in rows() if r["recipe_seed"] == 201 and (a is None or (r["a"], r["n"]) == (a, n))]


def size(row: dict) -> tuple[int, int]:
    return row["a"], row["n"]


def manifest() -> None:
    binaries = {"v0": ("IC_V0", V0_COMMIT), "base": ("IC_BASE", os.environ.get("IC_BASE_COMMIT")),
                "candidate": ("IC_CAND", os.environ.get("IC_CAND_COMMIT"))}
    doc = bench.host_manifest(RUNS / "host.json", {
        name: {"path_basename": arm(env).name, "sha256": bench.sha256(arm(env)), "built_from": commit}
        for name, (env, commit) in binaries.items()})
    r01 = json.loads((r01_runs() / "host.json").read_text())
    same = all(doc[k] == r01[k] for k in ("cpu_model", "cpu_flags_relevant", "logical_cores", "memory",
                                          "transparent_hugepage", "os"))
    (RUNS / "aa-source.json").write_text(json.dumps(
        {"host_matches_r01": same, "aa": "R01's" if same else "R02's own (step `aa`)"}, indent=1) + "\n")
    print(json.dumps({"host_matches_r01": same}, indent=1))


def per_summand(rep: dict) -> float:
    return stats.setup_phase_ns(rep)["collect"] / rep["counts"]["pass"]["collect"]["summands_scanned"]


def control() -> dict:
    """KIC_SCAN_SIMD=0 against the default on the baseline arm (PROTOCOL.md, "A control first")."""
    d = RUNS / "control"
    bench.interleave({"simd": arm("IC_BASE"), "scalar": arm("IC_BASE")}, m1(*CONTROL), ROUNDS, d,
                     envs={"simd": bench.ENV, "scalar": SCALAR_ENV})
    ratios = []
    for r in m1(*CONTROL):
        for k in range(1, ROUNDS + 1):
            a = bench.load(bench.figure_path(d / "scalar" / r["id"] / f"r{k}.price.json"))
            b = bench.load(bench.figure_path(d / "simd" / r["id"] / f"r{k}.price.json"))
            if a.get("status") == b.get("status") == "complete":
                ratios.append(per_summand(a) / per_summand(b))
    verdict = {"collect_ns_per_summand_scalar_over_simd": stats.geo_ci(ratios), "threshold": 1.3}
    verdict["confirmed"] = verdict["collect_ns_per_summand_scalar_over_simd"].get("geomean", 0) >= 1.3
    (d / "verdict.json").write_text(json.dumps(verdict, indent=1) + "\n")
    return verdict


def pin() -> dict:
    """Both arms' outputs on every suite row against v0's from R01 (amendment 3), untimed."""
    out = RUNS / "pin" / "pin.json"
    if out.exists():
        return json.loads(out.read_text())
    sys.path.insert(0, str(ROOT / "scripts"))
    import curve_id  # noqa: E402

    r01 = r01_runs()
    result = []
    for row in bench.slug_rows(bench.suite_rows("S") + bench.suite_rows("smoke")):
        if row["tier"] == "S":
            old = bench.load(bench.figure_path(r01 / "profile" / "v0" / row["suite_id"] / "r1.price.json"))
        else:
            old = bench.load(r01 / "smoke" / f"{row['suite_id']}.price.json")
        slug = bench.curve_slug(row["a"], row["n"])
        v0_name = curve_id.resolve(old.get("curve") or "")
        entry = {"row": row["id"], "suite_id": row["suite_id"], "slug": slug,
                 "v0_curve": old.get("curve"), "v0_curve_resolves_to_slug": (v0_name or {}).get("slug") == slug}
        for name, env in (("base", "IC_BASE"), ("candidate", "IC_CAND")):
            new = bench.untimed(arm(env), row, RUNS / "pin" / name / row["id"] / "pin.price.json")
            a, b = bench.outputs(new), bench.outputs(old)
            entry[name] = {"equal": a == b, "differs_in": sorted(k for k in a if a[k] != b[k]),
                           "curve": new.get("curve"), "curve_is_slug": new.get("curve") == slug}
        result.append(entry)
    arms = ("base", "candidate")
    doc = {"rows": result,
           "held": all(e[a]["equal"] for e in result for a in arms),
           "names_agree": all(e["v0_curve_resolves_to_slug"] and all(e[a]["curve_is_slug"] for a in arms)
                              for e in result)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def holdout_rows() -> list[dict]:
    """Seed 205 and targets T101/T102 per size, by suite v1's own construction."""
    spec = importlib.util.spec_from_file_location("suite_v1", bench.SUITE / "make_suite.py")
    suite = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(suite)
    out = []
    for (a, n), (columns, m, _) in suite.RECIPES.items():
        slug = bench.curve_slug(a, n)
        for i in HOLDOUT_TARGETS:
            text = suite.text(suite.params(a, n, columns, m, HOLDOUT_SEED, i))
            path = HERE / "holdouts" / slug / f"M5-T{i}.json"
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
            elif path.read_text() != text:
                raise SystemExit(f"{path} differs from its construction")
            out.append({"id": f"{slug}/M5-T{i}", "tier": "holdout", "a": a, "n": n, "r": suite.s20.R[(a, n)],
                        "recipe_seed": HOLDOUT_SEED, "target": i, "rho_seed": suite.RHO_SEED_BASE + i,
                        # Absolute, so bench.price_cmd's join with the suite leaves it as it is.
                        "params": str(path.resolve())})
    return out


def arms() -> dict[str, Path]:
    return {"base": arm("IC_BASE"), "cand": arm("IC_CAND")}


def compare() -> None:
    bench.interleave(arms(), rows(), ROUNDS, RUNS / "compare")


def cold_ratios(d: Path, rs: list[dict], rounds: int) -> list[float]:
    out = []
    for r in rs:
        for k in range(1, rounds + 1):
            a = bench.load(bench.figure_path(d / "base" / r["id"] / f"r{k}.price.json"))
            b = bench.load(bench.figure_path(d / "cand" / r["id"] / f"r{k}.price.json"))
            if a.get("status") == b.get("status") == "complete":
                out.append(stats.ic_cold_ns(a) / stats.ic_cold_ns(b))
    return out


def extend() -> list[str]:
    """Amendment 1: a size whose interval's half-width exceeds 5% gets five more rounds."""
    wide = []
    for key in sorted({size(r) for r in rows()}):
        rs = [r for r in rows() if size(r) == key]
        ci = stats.geo_ci(cold_ratios(RUNS / "compare", rs, ROUNDS))
        if "hi" in ci and ci["hi"] / ci["geomean"] - 1 > HALF_WIDTH_LIMIT:
            wide.append(rs[0]["id"].split("/")[0])
            bench.interleave(arms(), rs, EXTENDED_ROUNDS, RUNS / "compare")
    (RUNS / "compare" / "extended.json").write_text(json.dumps({"extended_sizes": wide}, indent=1) + "\n")
    return wide


def holdout() -> None:
    bench.interleave(arms(), holdout_rows(), ROUNDS, RUNS / "holdout")


def v0check() -> None:
    """v0 against v0' on M1's rows, five rounds: accounting, it gates nothing (amendment 3)."""
    bench.interleave({"v0": arm("IC_V0"), "base": arm("IC_BASE")}, m1(), ROUNDS, RUNS / "v0check")


def aa() -> None:
    """R02's own A/A on M1's rows, only when the host differs from R01's."""
    copy = arm("IC_BASE_COPY")
    if bench.sha256(copy) != bench.sha256(arm("IC_BASE")):
        raise SystemExit("IC_BASE_COPY is not a byte-identical copy of IC_BASE")
    bench.interleave({"A": arm("IC_BASE"), "A2": copy}, m1(), ROUNDS, RUNS / "aa")


def callgrind() -> None:
    """Amendment 2: the pipeline alone, `ic workflow`, under callgrind, both arms."""
    d = RUNS / "callgrind"
    d.mkdir(parents=True, exist_ok=True)
    for name, binary in arms().items():
        for a, n in CALLGRIND:
            row = next(r for r in m1(a, n) if r["target"] == 1)
            tag = f"{name}-{row['id'].replace('/', '-')}"
            stem = f"{tag}.callgrind.out"
            if (d / f"{tag}.phases.json").exists():
                continue
            work = Path(tempfile.mkdtemp(prefix="r02-workflow-"))
            cmd = ["valgrind", "--tool=callgrind", "--cache-sim=yes", "--I1=32768,8,64", "--D1=32768,8,64",
                   "--LL=2097152,16,64", f"--callgrind-out-file={d / stem}",
                   str(binary), "--json", "workflow", "--params", str(bench.SUITE / row["params"]),
                   "--dir", str(work)]
            with open(d / f"{tag}.workflow.json", "w") as out, open(d / f"{tag}.valgrind.log", "w") as log:
                subprocess.run(["taskset", "-c", bench.CPUS, *cmd], env=bench.ENV, stdout=out, stderr=log,
                               check=False)
            shutil.rmtree(work, ignore_errors=True)
            with open(d / f"{tag}.phases.json", "w") as f:
                subprocess.run([sys.executable, str(HERE.parents[1] / "harness" / "callgrind_phases.py"), str(d),
                                stem, "--top", "40"], stdout=f, check=False)


def rule() -> None:
    """§23's protocol at its six sizes with 64 targets, on the candidate."""
    s23 = ROOT / "research" / "ic_single_target_20260930"
    env = {**os.environ, "IC": str(arm("IC_CAND")), "IC_COMMIT": os.environ.get("IC_CAND_COMMIT", ""),
           "IC_RUNS": str(RUNS / "rule")}
    for a, n in [(1, 47), (0, 57), (0, 41), (0, 53), (1, 59), (0, 61)]:
        subprocess.run([sys.executable, str(s23 / "run.py"), "size", str(a), str(n)], env=env, check=True)


STEPS = {"manifest": manifest, "control": control, "pin": pin, "holdouts": holdout_rows, "compare": compare,
         "extend": extend, "holdout": holdout, "v0check": v0check, "callgrind": callgrind, "aa": aa,
         "rule": rule}


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
        raise SystemExit("an output differs from v0's; R02 stops (PROTOCOL.md, amendment 3)")
    if not json.loads((RUNS / "aa-source.json").read_text())["host_matches_r01"]:
        aa()
    compare()
    extend()
    holdout()
    v0check()
    callgrind()


if __name__ == "__main__":
    main()
