#!/usr/bin/env python3
"""R02b's steps, in the order PROTOCOL.md declares them.

    IC_BASE=<newest baseline> IC_BASE_COMMIT=<commit> \\
        IC_CAND=<candidate> IC_CAND_COMMIT=<commit> python3 run.py all
    python3 run.py <step>

The steps are manifest, pin, holdouts, compare, holdout, extend,
callgrind and rule. `all` runs the first seven in order and stops if the
pin fails; `rule` runs only when the round is accepted. The tests run before `all`, on the
candidate's tree, and are recorded beside the analysis.

Outputs go to runs/ (or IC_RUNS), one directory per row, named
`<slug>/<recipe>-<target>` (AGENTS.md §11). Nothing is overwritten: every
step resumes where it stopped.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402
import stats  # noqa: E402

ROOT = bench.ROOT
RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
R01_DIR = HERE.parent / "R01-baseline-v0"
ROUNDS, EXTENDED_ROUNDS, HALF_WIDTH_LIMIT = 5, 10, 0.03
TARGETS = [(1, 59), (0, 61)]
CALLGRIND = [(0, 61), (1, 59), (0, 41)]
HOLDOUTS = [(206, 103), (206, 104), (207, 105), (207, 106), (208, 107), (208, 108), (209, 109), (209, 110)]


def arm(name: str) -> Path:
    path = os.environ.get(name)
    if not path or not Path(path).exists():
        raise SystemExit(f"set {name}")
    return Path(path)


def arms() -> dict[str, Path]:
    return {"base": arm("IC_BASE"), "cand": arm("IC_CAND")}


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


def suite_rows() -> list[dict]:
    """`M1`'s 22 rows, and every suite row at the two target sizes: 34."""
    return [r for r in bench.slug_rows(bench.suite_rows("S"))
            if r["recipe_seed"] == 201 or (r["a"], r["n"]) in TARGETS]


def holdout_rows() -> list[dict]:
    """Eight fresh rows at each target size, by suite v1's own construction."""
    spec = importlib.util.spec_from_file_location("suite_v1", bench.SUITE / "make_suite.py")
    suite = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(suite)
    out = []
    for (a, n), (columns, m, _) in suite.RECIPES.items():
        if (a, n) not in TARGETS:
            continue
        slug = bench.curve_slug(a, n)
        for seed, i in HOLDOUTS:
            text = suite.text(suite.params(a, n, columns, m, seed, i))
            path = HERE / "holdouts" / slug / f"M{seed - 200}-T{i}.json"
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
            elif path.read_text() != text:
                raise SystemExit(f"{path} differs from its construction")
            out.append({"id": f"{slug}/M{seed - 200}-T{i}", "tier": "holdout", "a": a, "n": n,
                        "r": suite.s20.R[(a, n)], "recipe_seed": seed, "target": i,
                        "rho_seed": suite.RHO_SEED_BASE + i,
                        # Absolute, so bench.price_cmd's join with the suite leaves it as it is.
                        "params": str(path.resolve())})
    return out


def manifest() -> dict:
    binaries = {"base": ("IC_BASE", os.environ.get("IC_BASE_COMMIT")),
                "candidate": ("IC_CAND", os.environ.get("IC_CAND_COMMIT"))}
    doc = bench.host_manifest(RUNS / "host.json", {
        name: {"path_basename": arm(env).name, "sha256": bench.sha256(arm(env)), "built_from": commit}
        for name, (env, commit) in binaries.items()})
    r01 = json.loads((r01_runs() / "host.json").read_text())
    same = all(doc[k] == r01[k] for k in ("cpu_model", "cpu_flags_relevant", "logical_cores", "memory",
                                          "transparent_hugepage", "os"))
    out = {"host_matches_r01": same, "aa": "R01's" if same else "the round's own (not run by this script)"}
    (RUNS / "aa-source.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def pin() -> dict:
    """The candidate's outputs on every suite row against v0's from R01, untimed."""
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
        new = bench.untimed(arm("IC_CAND"), row, RUNS / "pin" / "candidate" / row["id"] / "pin.price.json")
        a, b = bench.outputs(new), bench.outputs(old)
        result.append({"row": row["id"], "suite_id": row["suite_id"], "slug": slug,
                       "equal": a == b, "differs_in": sorted(k for k in a if a[k] != b[k]),
                       "names_agree": (v0_name or {}).get("slug") == slug and new.get("curve") == slug})
    doc = {"rows": result, "held": all(e["equal"] for e in result),
           "names_agree": all(e["names_agree"] for e in result)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def compare() -> None:
    bench.interleave(arms(), suite_rows(), ROUNDS, RUNS / "compare")


def holdout() -> None:
    bench.interleave(arms(), holdout_rows(), ROUNDS, RUNS / "holdout")


def cold_ratios(d: Path, rs: list[dict], rounds: int) -> list[float]:
    out = []
    for r in rs:
        for k in range(1, rounds + 1):
            a = bench.load(bench.figure_path(d / "base" / r["id"] / f"r{k}.price.json"))
            b = bench.load(bench.figure_path(d / "cand" / r["id"] / f"r{k}.price.json"))
            if a.get("status") == b.get("status") == "complete":
                out.append(stats.ic_cold_ns(a) / stats.ic_cold_ns(b))
    return out


def extend() -> dict:
    """Rounds 6–10 for a set whose half-width at a target size exceeds 3%
    after five rounds. The test reads the interval's width only."""
    record = RUNS / "extended.json"
    if record.exists():
        done = json.loads(record.read_text())
    else:
        done = {}
        for name, d, rows in (("suite", RUNS / "compare", suite_rows()),
                              ("holdouts", RUNS / "holdout", holdout_rows())):
            for a, n in TARGETS:
                rs = [r for r in rows if (r["a"], r["n"]) == (a, n)]
                ci = stats.geo_ci(cold_ratios(d, rs, ROUNDS))
                if "hi" in ci and ci["hi"] / ci["geomean"] - 1 > HALF_WIDTH_LIMIT:
                    done.setdefault(name, []).append(bench.curve_slug(a, n))
        record.write_text(json.dumps(done, indent=1) + "\n")
    for name, d, rows in (("suite", RUNS / "compare", suite_rows()),
                          ("holdouts", RUNS / "holdout", holdout_rows())):
        rs = [r for r in rows if bench.curve_slug(r["a"], r["n"]) in done.get(name, [])]
        if rs:
            bench.interleave(arms(), rs, EXTENDED_ROUNDS, d)
    return done


def callgrind() -> None:
    """R02's control: `ic workflow` under callgrind, both arms, on `M1`'s
    first target at the two target sizes and at one narrow-kernel size."""
    d = RUNS / "callgrind"
    d.mkdir(parents=True, exist_ok=True)
    rows = suite_rows()
    for name, binary in arms().items():
        for a, n in CALLGRIND:
            row = next(r for r in rows if (r["a"], r["n"]) == (a, n) and r["recipe_seed"] == 201
                       and r["target"] == 1)
            tag = f"{name}-{row['id'].replace('/', '-')}"
            stem = f"{tag}.callgrind.out"
            if (d / f"{tag}.phases.json").exists():
                continue
            work = Path(tempfile.mkdtemp(prefix="r02b-workflow-"))
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
                                stem, "--top", "60"], stdout=f, check=False)


def rule() -> None:
    """§23's protocol at its six sizes with 64 targets, on the candidate."""
    s23 = ROOT / "research" / "ic_single_target_20260930"
    env = {**os.environ, "IC": str(arm("IC_CAND")), "IC_COMMIT": os.environ.get("IC_CAND_COMMIT", ""),
           "IC_RUNS": str(RUNS / "rule")}
    for a, n in [(1, 47), (0, 57), (0, 41), (0, 53), (1, 59), (0, 61)]:
        subprocess.run([sys.executable, str(s23 / "run.py"), "size", str(a), str(n)], env=env, check=True)


STEPS = {"manifest": manifest, "pin": pin, "holdouts": holdout_rows, "compare": compare, "holdout": holdout,
         "extend": extend, "callgrind": callgrind, "rule": rule}


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
    print(json.dumps(manifest(), indent=1), flush=True)
    doc = pin()
    if not (doc["held"] and doc["names_agree"]):
        raise SystemExit("an output differs from v0's, or a name disagrees; R02b stops (PROTOCOL.md)")
    holdout_rows()
    compare()
    holdout()
    print(json.dumps({"extended": extend()}, indent=1), flush=True)
    callgrind()


if __name__ == "__main__":
    main()
