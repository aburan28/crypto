#!/usr/bin/env python3
"""B4's measurement 5, read (PROTOCOL.md): every run's logarithm, checked
in the run and replayed here, `S` against rho's on the same point, and
every phase's cost and counts.

    IC_RUNS=<run tree> python3 analyse.py > f0.json

The replay is `[d]G = Q` in B1's generator arithmetic
(`../../conformance/v2/make_cases.py`), which shares nothing with the
tool, for both arms' certificates.  A run passes when its report is
`complete`, both arms' scalars replay, they agree, and a known
logarithm equals them.  The measurement passes when every run does.
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

_spec = importlib.util.spec_from_file_location(
    "conformance_v2_cases", HERE.parents[1] / "conformance" / "v2" / "make_cases.py")
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)


def as_int(v) -> int:
    return int(v, 16) if isinstance(v, str) and v.startswith("0x") else int(v)


def replay(doc: dict, scalar: int | None) -> dict:
    """`[scalar]G = Q` on the document's curve, and the known logarithm."""
    field, sub, target = doc["field"], doc["subgroup"], doc["target"]
    n = field["degree"]
    curve = v2.Curve(n, as_int(field["modulus"]), doc["curve"]["a"], 1)
    g = (as_int(sub["generator"]["x"]), as_int(sub["generator"]["y"]))
    known = as_int(target["known_log"]) if "known_log" in target else None
    q = (curve.mul(known, g) if known is not None
         else (as_int(target["point"]["x"]), as_int(target["point"]["y"])))
    replays = scalar is not None and curve.mul(scalar, g) == q
    return {"replays": replays, "known_log": known,
            "known_matches": None if known is None else scalar == known % as_int(sub["order"])}


def row(run_id: str, doc_path: Path, tree: Path) -> dict:
    out = tree / f"{run_id}.price.json"
    attempts = [run.attempt_path(out, k) for k in range(run.RETRIES + 1)]
    figure = next((a for a in attempts if a.exists() and run.bench.clean(a)
                   and run.bench.load(a).get("status") == "complete"), out)
    rep = run.bench.load(figure)
    doc = json.loads(doc_path.read_text())
    certs = rep.get("certificates") or {}
    ic_scalar = (certs.get("ic") or {}).get("scalar")
    rho_scalar = (certs.get("rho") or {}).get("scalar")
    ic_replay = replay(doc, int(ic_scalar) if ic_scalar is not None else None)
    rho_replay = replay(doc, int(rho_scalar) if rho_scalar is not None else None)
    median = rep.get("median") or {}
    reps = rep.get("repetitions") or [{}]
    record = run.bench.run_record(figure) or {}
    passed = (rep.get("status") == "complete" and ic_replay["replays"] and rho_replay["replays"]
              and ic_scalar == rho_scalar and ic_replay["known_matches"] is not False)
    return {
        "run": run_id,
        "document": str(doc_path.relative_to(HERE.parents[1])),
        "figure": figure.name,
        "attempts": sum(a.exists() for a in attempts),
        "exit_status": record.get("exit_status"),
        "contended": record.get("contended"),
        "status": rep.get("status"),
        "words": (rep.get("ic") or {}).get("words"),
        "scalar": ic_scalar,
        "verified_in_run": (rep.get("result") or {}).get("verified"),
        "ic_and_rho_agree": rep.get("ic_and_rho_agree"),
        "replay_ic": ic_replay,
        "replay_rho": rho_replay,
        "s_ic_cold": median.get("s_ic_cold"),
        "s_rho_cold": median.get("s_rho_cold"),
        "cold_ratio_ic_over_rho": median.get("cold_ratio_ic_over_rho"),
        "online_speedup": median.get("online_speedup"),
        "setup_phases_ns": reps[0].get("setup_phases_ns"),
        "ic_online_phases_ns": median.get("ic_online_phases_ns"),
        "counts": rep.get("counts"),
        "rho_counts": rep.get("rho_counts"),
        "pass": passed,
    }


def main() -> None:
    tree = Path(os.environ.get("IC_RUNS") or sys.exit("set IC_RUNS")) / "f0"
    rows = [row(run_id, doc, tree) for run_id, doc in run.documents()]
    print(json.dumps({"measurement": "B4's measurement 5: F0 at three words",
                      "runs": len(rows), "passed": sum(r["pass"] for r in rows),
                      "pass": all(r["pass"] for r in rows), "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
