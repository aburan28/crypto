#!/usr/bin/env python3
"""R04's steps, in the order PROTOCOL.md declares them.

    IC_DEFAULT=<binary> IC_PROBES=<binary> IC_COMMIT=<commit> python3 run.py all
    python3 run.py <step>   # manifest | pin | compare | aa

Both binaries are built from IC_COMMIT, the probes one with
`--features scan-probes`. `aa` runs only when the host differs from R01's,
and needs IC_DEFAULT_COPY, a byte-identical copy of IC_DEFAULT. Outputs
go to runs/ (or IC_RUNS), one directory per row, named
`<slug>/<recipe>-<target>` (AGENTS.md §11). Nothing is overwritten:
every step resumes where it stopped.
"""
from __future__ import annotations

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
ROUNDS = 3


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
    """`M1`'s 22 rows: two targets at each of the eleven sizes."""
    return [r for r in bench.slug_rows(bench.suite_rows("S")) if r["recipe_seed"] == 201]


def manifest() -> None:
    commit = os.environ.get("IC_COMMIT")
    doc = bench.host_manifest(RUNS / "host.json", {
        name: {"path_basename": arm(env).name, "sha256": bench.sha256(arm(env)), "built_from": commit,
               "features": features}
        for name, (env, features) in {"default": ("IC_DEFAULT", []),
                                      "probes": ("IC_PROBES", ["scan-probes"])}.items()})
    r01 = json.loads((r01_runs() / "host.json").read_text())
    same = all(doc[k] == r01[k] for k in ("cpu_model", "cpu_flags_relevant", "logical_cores", "memory",
                                          "transparent_hugepage", "os"))
    (RUNS / "aa-source.json").write_text(json.dumps(
        {"host_matches_r01": same, "aa": "R01's" if same else "R04's own (step `aa`)"}, indent=1) + "\n")
    print(json.dumps({"host_matches_r01": same}, indent=1))


def pin() -> dict:
    """Both arms' outputs on every row, untimed: they must be equal."""
    out = RUNS / "pin" / "pin.json"
    if out.exists():
        return json.loads(out.read_text())
    result = []
    for row in rows():
        a = bench.untimed(arm("IC_DEFAULT"), row, RUNS / "pin" / "default" / row["id"] / "pin.price.json")
        b = bench.untimed(arm("IC_PROBES"), row, RUNS / "pin" / "probes" / row["id"] / "pin.price.json")
        oa, ob = bench.outputs(a), bench.outputs(b)
        result.append({"row": row["id"], "equal": oa == ob, "differs_in": sorted(k for k in oa if oa[k] != ob[k]),
                       "probes_reported": "scan_probes" in b, "default_silent": "scan_probes" not in a})
    doc = {"rows": result,
           "held": all(r["equal"] for r in result),
           "probes_reported": all(r["probes_reported"] for r in result),
           "default_silent": all(r["default_silent"] for r in result)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc


def compare() -> None:
    bench.interleave({"default": arm("IC_DEFAULT"), "probes": arm("IC_PROBES")}, rows(), ROUNDS,
                     RUNS / "compare")


def aa() -> None:
    """R04's own A/A on `M1`'s rows, only when the host differs from R01's."""
    copy = arm("IC_DEFAULT_COPY")
    if bench.sha256(copy) != bench.sha256(arm("IC_DEFAULT")):
        raise SystemExit("IC_DEFAULT_COPY is not a byte-identical copy of IC_DEFAULT")
    bench.interleave({"A": arm("IC_DEFAULT"), "A2": copy}, rows(), 5, RUNS / "aa")


STEPS = {"manifest": manifest, "pin": pin, "compare": compare, "aa": aa}


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
    doc = pin()
    if not (doc["held"] and doc["probes_reported"] and doc["default_silent"]):
        raise SystemExit("the arms' outputs differ, or the probes are missing or leak; R04 stops (PROTOCOL.md)")
    if not json.loads((RUNS / "aa-source.json").read_text())["host_matches_r01"]:
        aa()
    compare()


if __name__ == "__main__":
    main()
