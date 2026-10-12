#!/usr/bin/env python3
"""R03's steps, in the order PROTOCOL.md declares them.

    IC_BASE=<base> IC_BASE_COMMIT=<commit> IC_CAND=<candidate> IC_CAND_COMMIT=<commit> \\
        python3 run.py all
    python3 run.py <step>   # manifest | pin | compare | holdout | aa

`aa` runs only when the host differs from R01's, and needs IC_BASE_COPY,
a byte-identical copy of IC_BASE.  Outputs go to runs/ (or IC_RUNS), one
directory per row, named `<slug>/<recipe>-<target>` (AGENTS.md §11).
Nothing is overwritten: every step resumes where it stopped.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()
R01_DIR = HERE.parent / "R01-baseline-v0"
R02_HOLDOUTS = HERE.parent / "R02-wide-tail-kernel" / "holdouts"
ROUNDS = 5
COMPOSITE = [(1, 45), (0, 57)]
HOLDOUT_TARGETS = (101, 102)


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


def size(row: dict) -> tuple[int, int]:
    return row["a"], row["n"]


def rows() -> list[dict]:
    """Every row at the two composite sizes, and `M1`'s rows elsewhere."""
    return [r for r in bench.slug_rows(bench.suite_rows("S"))
            if size(r) in COMPOSITE or r["recipe_seed"] == 201]


def manifest() -> None:
    binaries = {"base": ("IC_BASE", os.environ.get("IC_BASE_COMMIT")),
                "candidate": ("IC_CAND", os.environ.get("IC_CAND_COMMIT"))}
    doc = bench.host_manifest(RUNS / "host.json", {
        name: {"path_basename": arm(env).name, "sha256": bench.sha256(arm(env)), "built_from": commit}
        for name, (env, commit) in binaries.items()})
    r01 = json.loads((r01_runs() / "host.json").read_text())
    same = all(doc[k] == r01[k] for k in ("cpu_model", "cpu_flags_relevant", "logical_cores", "memory",
                                          "transparent_hugepage", "os"))
    (RUNS / "aa-source.json").write_text(json.dumps(
        {"host_matches_r01": same, "aa": "R01's" if same else "R03's own (step `aa`)"}, indent=1) + "\n")
    print(json.dumps({"host_matches_r01": same}, indent=1))


def pin() -> dict:
    """The candidate's outputs on every suite row against v0's from R01, untimed."""
    out = RUNS / "pin" / "pin.json"
    if out.exists():
        return json.loads(out.read_text())
    r01 = r01_runs()
    result = []
    for row in bench.slug_rows(bench.suite_rows("S") + bench.suite_rows("smoke")):
        if row["tier"] == "S":
            old = bench.load(bench.figure_path(r01 / "profile" / "v0" / row["suite_id"] / "r1.price.json"))
        else:
            old = bench.load(r01 / "smoke" / f"{row['suite_id']}.price.json")
        new = bench.untimed(arm("IC_CAND"), row, RUNS / "pin" / "candidate" / row["id"] / "pin.price.json")
        a, b = bench.outputs(new), bench.outputs(old)
        result.append({"row": row["id"], "suite_id": row["suite_id"], "equal": a == b,
                       "differs_in": sorted(k for k in a if a[k] != b[k]),
                       "curve_is_slug": new.get("curve") == bench.curve_slug(row["a"], row["n"])})
    doc = {"rows": result, "held": all(r["equal"] for r in result),
           "names_agree": all(r["curve_is_slug"] for r in result)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def holdout_rows() -> list[dict]:
    """R02's frozen holdout files at the two composite sizes."""
    spec = importlib.util.spec_from_file_location("suite_v1", bench.SUITE / "make_suite.py")
    suite = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(suite)
    out = []
    for a, n in COMPOSITE:
        slug = bench.curve_slug(a, n)
        for i in HOLDOUT_TARGETS:
            path = R02_HOLDOUTS / slug / f"M5-T{i}.json"
            if not path.exists():
                raise SystemExit(f"{path} is missing: R02's holdouts are this round's")
            out.append({"id": f"{slug}/M5-T{i}", "tier": "holdout", "a": a, "n": n, "r": suite.s20.R[(a, n)],
                        "recipe_seed": 205, "target": i, "rho_seed": suite.RHO_SEED_BASE + i,
                        # Absolute, so bench.price_cmd's join with the suite leaves it as it is.
                        "params": str(path.resolve())})
    return out


def arms() -> dict[str, Path]:
    return {"base": arm("IC_BASE"), "cand": arm("IC_CAND")}


def compare() -> None:
    bench.interleave(arms(), rows(), ROUNDS, RUNS / "compare")


def holdout() -> None:
    bench.interleave(arms(), holdout_rows(), ROUNDS, RUNS / "holdout")


def aa() -> None:
    """R03's own A/A on `M1`'s rows, only when the host differs from R01's."""
    copy = arm("IC_BASE_COPY")
    if bench.sha256(copy) != bench.sha256(arm("IC_BASE")):
        raise SystemExit("IC_BASE_COPY is not a byte-identical copy of IC_BASE")
    m1 = [r for r in bench.slug_rows(bench.suite_rows("S")) if r["recipe_seed"] == 201]
    bench.interleave({"A": arm("IC_BASE"), "A2": copy}, m1, ROUNDS, RUNS / "aa")


STEPS = {"manifest": manifest, "pin": pin, "holdouts": holdout_rows, "compare": compare, "holdout": holdout,
         "aa": aa}


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
    if not pin()["held"]:
        raise SystemExit("an output differs from v0's; R03 stops (PROTOCOL.md)")
    if not json.loads((RUNS / "aa-source.json").read_text())["host_matches_r01"]:
        aa()
    compare()
    holdout()


if __name__ == "__main__":
    main()
