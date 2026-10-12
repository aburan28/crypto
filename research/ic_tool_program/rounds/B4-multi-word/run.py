#!/usr/bin/env python3
"""B4's measurement 5 (PROTOCOL.md): F0 at three words.

    IC_B4=<ic binary> IC_B4_COMMIT=<sha> IC_RUNS=<run tree> python3 run.py manifest
    IC_B4=<ic binary> IC_B4_COMMIT=<sha> IC_RUNS=<run tree> python3 run.py f0

Each of design §4's six instances runs as `ic price` at F0 on two
targets, from B4's frozen documents in `../../conformance/v2-b4/params/`:
- the instance's case document (C088–C093), whose target is a known
  logarithm's multiple;
- its public point, `C0xx-T001`, whose logarithm the tool is not given.

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
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402

PARAMS = HERE.parents[1] / "conformance" / "v2-b4" / "params"
INSTANCES = (("C088", 127), ("C089", 137), ("C090", 151), ("C091", 157), ("C092", 173), ("C093", 179))
HOUR_S = 3600
RETRIES = 2


def documents() -> list[tuple[str, Path]]:
    """`(run id, document)`: each instance's two, keyed by its ICV1 slug."""
    out = []
    for cid, n in INSTANCES:
        known = PARAMS / f"{cid}-kic-three-word-n{n}.json"
        slug = re.search(r"icv1-[a-z0-9-]+", json.loads(known.read_text())["name"]).group(0)
        out.append((f"{slug}/{cid}-known", known))
        out.append((f"{slug}/{cid}-T001", PARAMS / f"{cid}-T001-n{n}.json"))
    return out


def binary() -> tuple[Path, str]:
    path, commit = os.environ.get("IC_B4"), os.environ.get("IC_B4_COMMIT")
    if not path or not Path(path).exists() or not commit:
        raise SystemExit("set IC_B4 and IC_B4_COMMIT")
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
        for k in range(RETRIES + 1):
            attempt = attempt_path(out, k)
            if not attempt.exists():
                attempt.parent.mkdir(parents=True, exist_ok=True)
                cmd = ["timeout", str(HOUR_S), str(ic), "price", "--params", str(doc), "--json",
                       "--out", str(attempt), "--repeats", "1", "--repeats-fast", "1"]
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
                              {"b4": {"path": str(ic), "built_from": commit, "sha256": bench.sha256(ic)}})
    print(json.dumps({k: doc[k] for k in ("cpu_model", "rustc", "binaries")}, indent=1))


if __name__ == "__main__":
    steps = {"manifest": manifest, "f0": f0}
    if len(sys.argv) != 2 or sys.argv[1] not in steps:
        raise SystemExit(f"usage: run.py {'|'.join(steps)}")
    steps[sys.argv[1]]()
