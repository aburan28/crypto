#!/usr/bin/env python3
"""Pre-outcome hash selection: field, one factor, and projected-column only."""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import signal
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925"
CORPUS = HERE.parent / "rotated_pdp_corpus_20260925"
sys.path.insert(0, str(PARENT))
sys.path.insert(0, str(CORPUS))
import gate  # noqa: E402
from producer import lifts  # noqa: E402 - #767 single-x lift only

DOMAIN = "ECC2K130-ROTATED-BETA-SWEEP-20260925-v1"
COUNTERS = 4096
CAP_SECONDS = 60
CAP_RSS = 128 * 1024 * 1024


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def orbit(field: gate.Field, beta: int) -> tuple[int, ...]:
    out, value = [], beta
    for _ in range(field.n):
        out.append(value)
        value = field.square(value)
    assert value == beta
    return tuple(out)


def column_counts(curve: gate.Curve, points: list[tuple[int, int]]) -> tuple[int, int]:
    projected = {curve.scalar(p, 4) for p in points}
    assert None in projected
    quotient = {None if p is None else min(p, curve.neg(p)) for p in projected}
    return len(projected), len(quotient)


def run_selection() -> dict:
    started, cpu_started = time.perf_counter(), time.process_time()
    field = gate.Field(19, [0, 1, 2, 5])
    field.rabin_prime_degree()
    curve = gate.Curve(field)
    baseline_orbit = set(orbit(field, 3))
    assert gate.rank(list(baseline_orbit)) == 19
    base = [3, orbit(field, 3)[6]]
    # The six-step conjugate is beta^(2^6), as in V_0 for m=6,d=2.
    assert base[1] == orbit(field, 3)[6]
    base_points = sorted(p for x in gate.all_x(base) for p in lifts(curve, x))
    assert len(base_points) == 7
    baseline_signed, baseline_quotient = column_counts(curve, base_points)
    attempts, primary = [], []
    seen_beta = set()
    selected_orbits = set(baseline_orbit)
    for counter in range(COUNTERS):
        payload = f"{DOMAIN}/beta/{counter}".encode("ascii")
        beta = 1 + int.from_bytes(hashlib.sha256(payload).digest(), "big") % ((1 << 19) - 1)
        row = {"counter": counter, "beta": beta}
        if beta in seen_beta:
            row["reason"] = "duplicate_beta"
        elif beta in selected_orbits:
            row["reason"] = "reference_or_selected_orbit"
        else:
            seen_beta.add(beta)
            conjugates = orbit(field, beta)
            row["normal_rank"] = gate.rank(list(conjugates))
            if row["normal_rank"] != 19:
                row["reason"] = "not_normal"
            else:
                row["trace"] = field.trace(beta)
                bases = gate.subspace_basis(list(conjugates), 6, 2)
                row["slice_ranks"] = [gate.rank(b) for b in bases]
                row["combined_rank"] = gate.rank([v for b in bases for v in b])
                assert row["trace"] == 1 and row["slice_ranks"] == [2] * 6
                assert row["combined_rank"] == 12
                points = sorted(p for x in gate.all_x(bases[0]) for p in lifts(curve, x))
                row["f0_size"] = len(points)
                if len(points) != 7:
                    row["reason"] = "factor_size"
                else:
                    signed, quotient = column_counts(curve, points)
                    row["projected_signed"] = signed
                    row["projected_sign_quotient"] = quotient
                    row["reason"] = "equal_column" if signed == baseline_signed else "unequal_column"
                    if signed == baseline_signed:
                        primary.append({"counter": counter, "beta": beta,
                                        "projected_signed": signed,
                                        "projected_sign_quotient": quotient,
                                        "selection": "primary"})
                        selected_orbits.update(conjugates)
        attempts.append(row)
        if len(primary) == 4:
            break
    selected = list(primary)
    fallback_used = len(primary) < 3
    if fallback_used:
        for row in attempts:
            if row.get("reason") != "unequal_column" or len(selected) == 4:
                continue
            beta = row["beta"]
            if beta in selected_orbits:
                continue
            selected.append({"counter": row["counter"], "beta": beta,
                             "projected_signed": row["projected_signed"],
                             "projected_sign_quotient": row["projected_sign_quotient"],
                             "selection": "fallback"})
            selected_orbits.update(orbit(field, beta))
    outcome = {"status": "success" if len(selected) >= 3 else "insufficient_candidates",
               "domain": DOMAIN, "scan_limit": COUNTERS,
               "reference_beta": 3,
               "reference_f0_size": len(base_points),
               "reference_projected_signed": baseline_signed,
               "reference_projected_sign_quotient": baseline_quotient,
               "fallback_used": fallback_used,
               "selected": selected, "attempts": attempts,
               "field_operations": dict(field.operations),
               "curve_operations": dict(curve.operations),
               "wall_seconds": time.perf_counter() - started,
               "cpu_seconds": time.process_time() - cpu_started,
               "peak_rss_bytes": rss()}
    assert outcome["wall_seconds"] <= CAP_SECONDS and outcome["peak_rss_bytes"] <= CAP_RSS
    return outcome


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    def expired(_signal, _frame):
        raise TimeoutError(f"preflight {CAP_SECONDS}s wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        result = run_selection()
    except Exception as error:
        result = {"status": "failed", "error": repr(error),
                  "wall_seconds": time.perf_counter() - started,
                  "peak_rss_bytes": rss()}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    if result["status"] == "failed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
