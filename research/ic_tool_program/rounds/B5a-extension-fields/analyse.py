#!/usr/bin/env python3
"""B5a's measurement 5, read (PROTOCOL.md): every run's logarithm, checked
in the run and replayed here, each arm's `S`, and every phase's cost and
counts.

    IC_RUNS=<run tree> python3 analyse.py > f0.json

The replay is `[d]G = Q` in the instance generator's arithmetic for
`GF(p^k)` (`instances.py`), which shares nothing with the tool, for every
arm's certificate.  A run passes when its report is `complete`, every
arm's scalar replays, the arms agree, and a known logarithm equals them.
The measurement passes when every run does.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run  # noqa: E402

_spec = importlib.util.spec_from_file_location("b5a_instances", HERE / "instances.py")
inst = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(inst)


def point(v: dict) -> tuple:
    return tuple(int(c) for c in v["x"]), tuple(int(c) for c in v["y"])


def replay(doc: dict, scalar: int | None) -> dict:
    """`[scalar]G = Q` on the document's curve, and the known logarithm."""
    field, sub, target = doc["field"], doc["subgroup"], doc["target"]
    F = inst.Fpk(int(field["p"]), [int(c) for c in field["modulus"]])
    E = inst.Curve(F, tuple(int(c) for c in doc["curve"]["a"]), tuple(int(c) for c in doc["curve"]["b"]))
    G = point(sub["generator"])
    known = int(target["known_log"]) if "known_log" in target else None
    Q = E.mul(known, G) if known is not None else point(target["point"])
    replays = scalar is not None and E.mul(scalar, G) == Q
    return {"replays": replays, "known_log": known,
            "known_matches": None if known is None else scalar == known % int(sub["order"])}


def row(run_id: str, doc_path: Path, tree: Path) -> dict:
    out = tree / f"{run_id}.price.json"
    attempts = [run.attempt_path(out, k) for k in range(run.RETRIES + 1)]
    figure = next((a for a in attempts if a.exists() and run.bench.clean(a)
                   and run.bench.load(a).get("status") == "complete"), out)
    rep = run.bench.load(figure)
    doc = json.loads(doc_path.read_text())
    certs = rep.get("certificates") or {}
    scalars = {arm: (certs.get(arm) or {}).get("scalar") for arm in ("ic", "rho")}
    scalars = {arm: s for arm, s in scalars.items() if s is not None}
    if not scalars and (rep.get("result") or {}).get("scalar") is not None:
        scalars = {"rho": rep["result"]["scalar"]}
    replays = {arm: replay(doc, int(s)) for arm, s in scalars.items()}
    median = rep.get("median") or {}
    record = run.bench.run_record(figure) or {}
    passed = (rep.get("status") == "complete" and bool(replays)
              and all(r["replays"] and r["known_matches"] is not False for r in replays.values())
              and len({str(s) for s in scalars.values()}) == 1)
    return {
        "run": run_id,
        "document": str(doc_path.relative_to(HERE.parents[1])),
        "figure": figure.name,
        "attempts": sum(a.exists() for a in attempts),
        "exit_status": record.get("exit_status"),
        "contended": record.get("contended"),
        "status": rep.get("status"),
        "route": rep.get("route"),
        "scalars": scalars,
        "verified_in_run": (rep.get("result") or {}).get("verified"),
        "replays": replays,
        "s_ic_cold": median.get("s_ic_cold"),
        "s_rho_cold": median.get("s_rho_cold"),
        "cold_ratio_ic_over_rho": median.get("cold_ratio_ic_over_rho"),
        "counts": rep.get("counts"),
        "rho_counts": rep.get("rho_counts"),
        "pass": passed,
    }


def main() -> None:
    tree = Path(os.environ.get("IC_RUNS") or sys.exit("set IC_RUNS")) / "f0"
    rows = [row(run_id, doc, tree) for run_id, doc in run.documents()]
    print(json.dumps({"measurement": "B5a's measurement 5: F0 on extension fields",
                      "runs": len(rows), "passed": sum(r["pass"] for r in rows),
                      "pass": all(r["pass"] for r in rows), "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
