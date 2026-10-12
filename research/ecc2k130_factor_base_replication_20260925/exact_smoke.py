#!/usr/bin/env python3
"""Conditional bounded ECC2K-130 degree-263 representation smoke, no DLP claim."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_factor_base_pilot_20260924"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_dual_transport_20260925"))
import run as pilot  # noqa: E402
from dual_transport import DualTransport  # noqa: E402
from fastfield import IRR131  # noqa: E402
from relations import (CHALLENGE_PX, CHALLENGE_PY, CHALLENGE_QX,
                       CHALLENGE_QY, CHALLENGE_ELL)  # noqa: E402

BASE_SIZE = 16
TARGETS = 32
DEGREE = 263
COFACTOR = 4
WALL_CAP = 300
RSS_CAP_BYTES = 2 * (1 << 30)
VARIANTS = ("original", "transported", "descendant_native", "pullback")


def h_int(label):
    return int.from_bytes(hashlib.sha256(label.encode()).digest(), "big")


def phase(meter, ledger, name, fn):
    before = meter.snapshot()
    result = fn()
    ledger[name] = meter.delta(before, meter.snapshot())
    return result


def actual_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value * 1024 if sys.platform.startswith("linux") else value


def check_cap(start):
    if time.monotonic() - start > WALL_CAP:
        raise TimeoutError("line exceeded 300 wall seconds")
    if actual_rss_bytes() > RSS_CAP_BYTES:
        raise MemoryError("line exceeded 2 GiB peak RSS")


def point_hex(point):
    return None if point is None else [hex(point[0]), hex(point[1])]


def from_hex(pair):
    return tuple(int(value, 16) for value in pair)


def prefix_base(E, label, start):
    base, seen, meta = [], set(), Counter()
    for x in range(4096):
        check_cap(start)
        meta["x_scanned"] += 1
        for P in E.points_over(x):
            meta["cofactor_projections"] += 1
            Q = E.mul(P, COFACTOR)
            if Q is None:
                meta["infinity"] += 1
            elif Q in seen:
                meta["duplicate"] += 1
            else:
                seen.add(Q)
                base.append(Q)
                if len(base) == BASE_SIZE:
                    return base, dict(meta)
    raise AssertionError(f"{label} only yielded {len(base)} useful points")


def independent_rank(rows, modulus):
    matrix = [r[:] for r in rows]
    cursor = 0
    for col in range(BASE_SIZE + 1):
        pivot = next((j for j in range(cursor, len(matrix)) if matrix[j][col] % modulus), None)
        if pivot is None:
            continue
        matrix[cursor], matrix[pivot] = matrix[pivot], matrix[cursor]
        inverse = pow(matrix[cursor][col], -1, modulus)
        matrix[cursor] = [z * inverse % modulus for z in matrix[cursor]]
        for j in range(cursor + 1, len(matrix)):
            multiplier = matrix[j][col]
            if multiplier:
                matrix[j] = [(a - multiplier * b) % modulus
                             for a, b in zip(matrix[j], matrix[cursor])]
        cursor += 1
        if cursor == len(matrix):
            break
    return cursor


def scan(E, base, targets, r):
    table = pilot.pair_table(E, base)
    cases, rows = [], []
    for u, v, T in targets:
        hits = table.get(T, [])
        independent = False
        if hits:
            a, b = hits[0]
            row = [0] * (BASE_SIZE + 1)
            row[a] += 1
            row[b] += 1
            row[-1] = -v % r
            before = independent_rank(rows, r)
            rows.append(row)
            independent = independent_rank(rows, r) > before
        cases.append({"witness_count": len(hits),
                      "first_witness": list(hits[0]) if hits else None,
                      "independent": independent})
    return {"hits": sum(c["witness_count"] > 0 for c in cases),
            "misses": sum(c["witness_count"] == 0 for c in cases),
            "independent_rank": independent_rank(rows, r),
            "distinct_pair_sums": len(table), "pair_entries": sum(map(len, table.values())),
            "cases": cases}


def single_line(item, basis_v, P, Q, field_modulus):
    started = time.monotonic()
    meter, ledger = pilot.Meter(), {}
    before = meter.snapshot()
    F = pilot.CountedField(meter, deg=131, irr=field_modulus)
    F.frobenius(1, 130)
    ledger["field_setup"] = meter.delta(before, meter.snapshot())
    E0, Twist = pilot.CountedCurve(meter, F, 0, 1), pilot.CountedCurve(meter, F, 1, 1)
    r = CHALLENGE_ELL
    assert E0.on_curve(P) and E0.on_curve(Q)
    assert E0.mul(P, r) is None and E0.mul(Q, r) is None
    check_cap(started)
    H = from_hex(item["kernel_generator_on_twist"])
    assert Twist.mul(H, DEGREE) is None
    forward = phase(meter, ledger, "forward_setup", lambda: pilot.BinaryVeluMap.from_generator(
        E0, Twist, H, DEGREE))
    D = phase(meter, ledger, "dual_setup", lambda: DualTransport(
        E0, Twist, H, basis_v, DEGREE, P))
    assert (forward.codomain.a, forward.codomain.b) == (D.codomain.a, D.codomain.b)
    assert D.compose(P) == E0.mul(P, DEGREE)
    assert D.compose(Q) == E0.mul(Q, DEGREE)
    E1 = pilot.CountedCurve(meter, F, forward.codomain.a, forward.codomain.b)
    assert E1.b == int(item["codomain_b"], 16)
    check_cap(started)
    B0, meta0 = phase(meter, ledger, "source_base",
                       lambda: prefix_base(E0, "source", started))
    B1, meta1 = phase(meter, ledger, "native_base",
                       lambda: prefix_base(E1, "native", started))
    Btransport = phase(meter, ledger, "base_transport",
                       lambda: [forward(S) for S in B0])
    inv_degree = pow(DEGREE, -1, r)
    Bpullback = phase(meter, ledger, "native_dual_pullback",
                       lambda: [E0.mul(D.dual(S), inv_degree) for S in B1])
    assert all(forward(src) == leaf for src, leaf in zip(Bpullback, B1))
    assert len(set(Btransport)) == len(set(Bpullback)) == BASE_SIZE
    # Native leaf orbit action is phi(tau(phi^-1(S))). It is a charged
    # representative neighbor per selected point, separate from pair scanning.
    orbit_neighbors = phase(meter, ledger, "native_leaf_orbit_action",
        lambda: [forward(E0.frobenius(S)) for S in Bpullback])
    assert all(E1.on_curve(S) and E1.mul(S, r) is None for S in orbit_neighbors)
    P1, Q1 = phase(meter, ledger, "codomain_challenge_transport",
                   lambda: (forward(P), forward(Q)))
    assert E1.mul(P1, r) is None and E1.mul(Q1, r) is None
    check_cap(started)
    label = f"{item['line'][0]},{item['line'][1]}"
    coefficients = [(h_int(f"degree263-smoke-v1|{label}|{i}|u") % r,
                     1 + h_int(f"degree263-smoke-v1|{label}|{i}|v") % (r - 1))
                    for i in range(TARGETS)]
    targets0, targets1, costs0, costs1 = [], [], [], []
    for u, v in coefficients:
        check_cap(started)
        before = meter.snapshot()
        targets0.append((u, v, E0.add(E0.mul(P, u), E0.mul(Q, v))))
        costs0.append(meter.delta(before, meter.snapshot()))
        before = meter.snapshot()
        targets1.append((u, v, E1.add(E1.mul(P1, u), E1.mul(Q1, v))))
        costs1.append(meter.delta(before, meter.snapshot()))
    phase(meter, ledger, "target_covariance_audit",
          lambda: [(_ for _ in ()).throw(AssertionError("target map mismatch"))
                   if forward(src[2]) != dst[2] else None
                   for src, dst in zip(targets0, targets1)])
    bases = {"original": B0, "transported": Btransport,
             "descendant_native": B1, "pullback": Bpullback}
    results, scan_costs = {}, {}
    for name in VARIANTS:
        E = E1 if name in ("transported", "descendant_native") else E0
        stream = targets1 if E is E1 else targets0
        results[name] = phase(meter, scan_costs, name,
                              lambda E=E, base=bases[name], stream=stream:
                              scan(E, base, stream, r))
    assert [c["witness_count"] for c in results["original"]["cases"]] == [
        c["witness_count"] for c in results["transported"]["cases"]]
    assert [c["witness_count"] for c in results["descendant_native"]["cases"]] == [
        c["witness_count"] for c in results["pullback"]["cases"]]
    return {"line": item["line"], "status": "PASS", "codomain_b": hex(E1.b),
            "raw_reverse_sign": D.raw_reverse_scalar_sign,
            "base_meta": {"source": meta0, "native": meta1},
            "bases": {name: [point_hex(S) for S in base] for name, base in bases.items()},
            "orbit_neighbor_count": len(orbit_neighbors),
            "coefficient_pairs": [[u, v] for u, v in coefficients],
            "results": results, "phase_costs": ledger, "scan_costs": scan_costs,
            "source_target_prefix_costs": costs0, "leaf_target_prefix_costs": costs1,
            "wall_seconds": time.monotonic() - started,
            "peak_rss_bytes": actual_rss_bytes()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy", type=Path, required=True)
    ap.add_argument("--replay", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    assert not args.out.exists()
    toy, replay = json.loads(args.toy.read_text()), json.loads(args.replay.read_text())
    assert replay["status"] == "PASS" and replay["input_sha256"] == hashlib.sha256(args.toy.read_bytes()).hexdigest()
    assert all(toy["controls"][s][h][key] for s in toy["controls"] for h in toy["controls"][s]
               for key in toy["controls"][s][h])
    assert all(toy["variants"][s][h][v]["rank"] == BASE_SIZE + 1
               for s in toy["variants"] for h in toy["variants"][s] for v in VARIANTS)
    preflight_path = ROOT / "research/ecc2k130_direction_review_20260924/twist_torsion_results.json"
    preflight = json.loads(preflight_path.read_text())
    run = next(x for x in preflight["runs"] if x["seed"] == 20260924)
    reps = {tuple(x["line"]): x for x in run["representative_x_maps"]}
    V = from_hex(run["basis_on_twist"][1])
    P, Q = (CHALLENGE_PX, CHALLENGE_PY), (CHALLENGE_QX, CHALLENGE_QY)
    output = {"schema": "ecc2k130-degree263-bounded-native-transport-smoke-v1",
              "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "status": "PASS", "platform": platform.platform(),
              "toy_sha256": hashlib.sha256(args.toy.read_bytes()).hexdigest(),
              "replay_sha256": hashlib.sha256(args.replay.read_bytes()).hexdigest(),
              "preflight_sha256": hashlib.sha256(preflight_path.read_bytes()).hexdigest(),
              "preflight_shared_wall_seconds_unmetered": run["elapsed_seconds"],
              "public_challenge_input": {
                  "repository_path": "research/ecc2k130_relations/relations.py",
                  "repository_sha256": hashlib.sha256((ROOT / "research/ecc2k130_relations/relations.py").read_bytes()).hexdigest(),
                  "primary_parameters_url": "https://www.certicom.com/en/curves-list",
                  "independent_paper_url": "https://www.ecc-challenge.info/anon.pdf",
                  "P": point_hex(P), "Q": point_hex(Q),
                  "subgroup_order": str(CHALLENGE_ELL),
                  "curve_and_subgroup_membership_checked_per_line": True},
              "caps": {"targets_per_line": TARGETS, "wall_seconds_per_line": WALL_CAP,
                       "rss_bytes_per_line": RSS_CAP_BYTES}, "lines": []}
    for line in ((1, 0), (1, 4)):
        try:
            result = single_line(reps[line], V, P, Q, IRR131)
        except (TimeoutError, MemoryError, AssertionError, ArithmeticError) as error:
            result = {"line": line, "status": "STOP", "reason": repr(error),
                      "peak_rss_bytes": actual_rss_bytes()}
            output["status"] = "STOP"
        output["lines"].append(result)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": output["status"],
                      "lines": [{"line": x["line"], "status": x["status"],
                                 "hits": {k: v["hits"] for k, v in x.get("results", {}).items()},
                                 "seconds": x.get("wall_seconds"),
                                 "reason": x.get("reason")}
                                for x in output["lines"]], "output": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()
