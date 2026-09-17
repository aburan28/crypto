#!/usr/bin/env python3
"""Cite frozen control/candidate receipts. Does not invent rates."""
from __future__ import annotations

import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLOOR_LO = 23.0
FLOOR_HI = 25.0
HISTORICAL_BENCH = 14.637530


def samples(obj: dict) -> list[float]:
    rows = obj.get("samples") or obj.get("benchmark", {}).get("samples") or []
    rates = []
    for row in rows:
        if not row.get("valid", False):
            raise SystemExit("receipt contains an invalid sample")
        rate = float(row["rate"])
        if rate > 1000:
            rate /= 1000.0
        rates.append(rate)
    if not rates and obj.get("rate"):
        rate = float(obj["rate"])
        rates.append(rate / 1000.0 if rate > 1000 else rate)
    if not rates:
        raise SystemExit("receipt has no rates")
    return rates


def main() -> None:
    control_path = HERE / "control.json"
    candidate_path = HERE / "candidate.json"
    if not control_path.exists() or not candidate_path.exists():
        print("receipts not written yet; nothing to summarize")
        return
    control = json.loads(control_path.read_text())
    candidate = json.loads(candidate_path.read_text())
    if not control.get("valid") or not candidate.get("valid"):
        raise SystemExit("one or both receipts are invalid")
    cr = samples(control)
    xr = samples(candidate)
    cmed = statistics.median(cr)
    xmed = statistics.median(xr)
    ratio = xmed / cmed
    hist = xmed / HISTORICAL_BENCH
    floor_lo = xmed / FLOOR_HI
    floor_hi = xmed / FLOOR_LO
    pairs = None
    if len(cr) == len(xr):
        pairs = all(a < b for a, b in zip(cr, xr))
    accept = ratio >= 1.01 and (pairs is True)
    summary = {
        "unit": "B scalar updates / s",
        "class": "engineering",
        "controlRatesB": cr,
        "candidateRatesB": xr,
        "controlMedianB": cmed,
        "candidateMedianB": xmed,
        "ratioToMatchedControl": ratio,
        "ratioToHistorical14_637530": hist,
        "ratioToOneAddFloor23": floor_hi,
        "ratioToOneAddFloor25": floor_lo,
        "everyPairFaster": pairs,
        "acceptance1pctAllPairs": accept,
        "promoteDefault": accept,
        "controlIdentity": {
            "packedTopClmad": control.get("packedTopClmad"),
            "expectedPackedTopClmad": control.get("expectedPackedTopClmad"),
        },
        "candidateIdentity": {
            "packedTopClmad": candidate.get("packedTopClmad"),
            "expectedPackedTopClmad": candidate.get("expectedPackedTopClmad"),
        },
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
