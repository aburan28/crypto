#!/usr/bin/env python3
"""B5a's measurement 5 (PROTOCOL.md): F0 on extension fields.

    IC_B5A=<ic binary> IC_B5A_COMMIT=<sha> IC_RUNS=<run tree> python3 run.py manifest
    IC_B5A=<ic binary> IC_B5A_COMMIT=<sha> IC_RUNS=<run tree> python3 run.py f0

Each instance of design §4 runs as `ic price` at F0 on two targets, from
B5a's frozen documents in `../../conformance/v2-b5a/params/`:
- `<id>-known.json`, whose target is a known logarithm's multiple;
- `<id>-T001.json`, a public point whose logarithm the tool is not given.

`G1`–`G3` are paired (`ic-gaudry-cubic` against `rho-negation`); the
others are `solve: rho`, as their documents say. C050's document, B2's,
runs on its own known-answer target.

Each run is one process through the programme's runner
(`../../harness/bench.py`'s `launch`: `tools/isolated_bench.py run
--wait --cpus 2`, `RAYON_NUM_THREADS=1`, after PSI has fallen), with an
hour (`timeout 3600`).  A contended or failed attempt is kept and the
run made again, at most twice; the first clean, complete attempt is the
run's figure.  Nothing is ever overwritten, so each command resumes
where the last stopped.  `analyse.py` reads the run tree.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402

PARAMS = HERE.parents[1] / "conformance" / "v2-b5a" / "params"
C050 = HERE.parents[1] / "conformance" / "v2-b2" / "params" / "C050-cubic-extension.json"
INSTANCES = ("G1", "G2", "G3", "E2", "E5", "E11", "B2")
HOUR_S = 3600
RETRIES = 2


def documents() -> list[tuple[str, Path]]:
    """`(run id, document)`: each instance's two, keyed by its ICV1 slug, then C050's."""
    recs = {r["id"]: r for r in json.loads((HERE / "instances.json").read_text())["instances"]}
    out = []
    for iid in INSTANCES:
        slug = recs[iid]["slug"]
        out.append((f"{slug}/{iid}-known", PARAMS / f"{iid}-known.json"))
        out.append((f"{slug}/{iid}-T001", PARAMS / f"{iid}-T001.json"))
    out.append(("icv1-fp10k3-t41822-dfa3991d/C050-known", C050))
    return out


def binary() -> tuple[Path, str]:
    path, commit = os.environ.get("IC_B5A"), os.environ.get("IC_B5A_COMMIT")
    if not path or not Path(path).exists() or not commit:
        raise SystemExit("set IC_B5A and IC_B5A_COMMIT")
    return Path(path), commit


def runs() -> Path:
    path = os.environ.get("IC_RUNS")
    if not path:
        raise SystemExit("set IC_RUNS")
    return Path(path)


def attempt_path(out: Path, k: int) -> Path:
    return out if k == 0 else out.with_name(bench.stem(out).name + f"-retry{k}.price.json")


def f0() -> None:
    ic, _ = binary()
    tree = runs() / "f0"
    for run_id, doc in documents():
        out = tree / f"{run_id}.price.json"
        paired = json.loads(doc.read_text())["method"]["solve"] == "paired"
        for k in range(RETRIES + 1):
            attempt = attempt_path(out, k)
            if not attempt.exists():
                attempt.parent.mkdir(parents=True, exist_ok=True)
                cmd = ["timeout", str(HOUR_S), str(ic), "price", "--params", str(doc), "--json",
                       "--out", str(attempt)] + (["--repeats", "1", "--repeats-fast", "1"] if paired else [])
                bench.launch(cmd, attempt, tree)
            rep = bench.load(attempt)
            if bench.clean(attempt) and rep.get("status") == "complete":
                break
        result = rep.get("result", {})
        print(f"{run_id}: {rep.get('status')} scalar {result.get('scalar')} "
              f"verified {result.get('verified')}", flush=True)


def manifest() -> None:
    ic, commit = binary()
    tree = runs()
    tree.mkdir(parents=True, exist_ok=True)
    doc = bench.host_manifest(tree / "host.json",
                              {"b5a": {"path": str(ic), "built_from": commit, "sha256": bench.sha256(ic)}})
    print(json.dumps({k: doc[k] for k in ("cpu_model", "rustc", "binaries")}, indent=1))


if __name__ == "__main__":
    steps = {"manifest": manifest, "f0": f0}
    if len(sys.argv) != 2 or sys.argv[1] not in steps:
        raise SystemExit(f"usage: run.py {'|'.join(steps)}")
    steps[sys.argv[1]]()
