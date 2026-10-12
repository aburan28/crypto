#!/usr/bin/env python3
"""Exact n19 phase census from the pinned complete projected histogram."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
ARCHIVE = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
GATE = NOTES / "rotated_subspace_support_20260925/gate.py"
POINTS = NOTES / "rotated_four_base_joint_rank_20260925/inputs/point_only.json"
INPUT = HERE / "INPUT.json"
TORSION = (None, (0, 1), (1, 0), (1, 1))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def point(value):
    return None if value is None else tuple(value)


def pjson(value):
    return None if value is None else list(value)


def arithmetic():
    spec = importlib.util.spec_from_file_location("phase_producer_arithmetic", GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tar_json(tar: tarfile.TarFile, name: str):
    member = tar.getmember(name)
    assert member.isfile() and 0 < member.size < 40_000_000
    stream = tar.extractfile(member)
    assert stream is not None
    return json.load(stream)


def tar_jsonl(tar: tarfile.TarFile, name: str):
    member = tar.getmember(name)
    assert member.isfile() and 0 < member.size < 40_000_000
    stream = tar.extractfile(member)
    assert stream is not None
    return [json.loads(line) for line in stream]


def histogram(curve, manifest: dict):
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        prefix = "raw/n19-m6/"
        factors = [[point(value) for value in slot]
                   for slot in tar_json(tar, prefix + "factors.json")]
        rows = tar_jsonl(tar, prefix + "projected_histogram.jsonl")
    assert len(factors) == manifest["m"]
    assert all(len(slot) == manifest["physical_points_per_slot"] for slot in factors)
    assert all(curve.on_curve(p) for slot in factors for p in slot)
    assert all(set(curve.tau(p) for p in factors[i]) == set(factors[i + 1])
               for i in range(5))
    result = {}
    total = 0
    for row in rows:
        target = point(row["point"])
        assert target not in result and row["count"] > 0
        assert len(row["witness_indices"]) == 6
        assert all(0 <= index < 7 for index in row["witness_indices"])
        result[target] = row
        total += row["count"]
    assert len(result) == manifest["reference_support"]
    assert total == manifest["tuple_count"]
    assert None in result
    return factors, result


def column_terms(curve, factors, q: int, lam: int):
    slots = []
    all_columns = set()
    for i, slot in enumerate(factors):
        terms = []
        for source in slot:
            back = source
            for _ in range((19 - i) % 19):
                back = curve.tau(back)
            assert back in factors[0]
            projected = curve.scalar(back, 4)
            if projected is None:
                assert back in TORSION and source in TORSION
                terms.append((None, 0))
                continue
            canonical = min(projected, curve.neg(projected))
            sign = 1 if canonical == projected else -1
            coeff = sign * pow(lam, i, q) % q
            assert curve.scalar(source, 4) == curve.scalar(canonical, coeff)
            terms.append((canonical, coeff))
            all_columns.add(canonical)
        slots.append(terms)
    base_columns = {column for column, _ in slots[0] if column is not None}
    assert len(base_columns) == 3 and all_columns == base_columns
    return slots, sorted(base_columns)


def exact_index(curve, h, q: int):
    index = {}
    current = None
    for j in range(q):
        assert current not in index
        index[current] = j
        current = curve.add(current, h)
    assert current is None and len(index) == q
    return index


def first_hits(support: bytearray, powers: list[int], order: list[int]) -> bytes:
    q = len(support)
    found = bytearray([255]) * q
    for j in range(q):
        for position, phase in enumerate(order):
            if support[(powers[phase] * j) % q]:
                found[j] = position
                break
    return bytes(found)


def cap_summary(first: bytes, order: list[int], caps: list[int], q: int):
    result = {}
    for cap in caps:
        assert cap <= len(order)
        histogram = [0] * cap
        probes = 0
        supported = 0
        for position in first:
            if position < cap:
                supported += 1
                histogram[position] += 1
                probes += position + 1
            else:
                probes += cap
        misses = q - supported
        # The second identity avoids any independence assumption.
        prefixes = [sum(1 for position in first if position < t) for t in range(cap)]
        assert probes == sum(q - covered for covered in prefixes)
        result[str(cap)] = {
            "support": supported,
            "nonzero_support": supported - int(first[0] < cap),
            "misses": misses,
            "first_hit_position_counts": histogram,
            "phase_at_position": order[:cap],
            "oracle_probes": probes,
            "probe_mean_denominator": q,
        }
    return result


def certificate(curve, qpoint, position: int, order: list[int],
                factors, rows, terms, q: int, lam: int):
    if position == 255:
        return {"status": "miss", "attempts": 19}
    phase = order[position]
    shifted = qpoint
    for _ in range(phase):
        shifted = curve.tau(shifted)
    assert shifted == curve.scalar(qpoint, pow(lam, phase, q))
    target_r = curve.scalar(shifted, 4)
    archive = rows[target_r]
    indices = archive["witness_indices"]
    total = None
    coefficients = {}
    for i, index in enumerate(indices):
        total = curve.add(total, factors[i][index])
        column, coeff = terms[i][index]
        if column is not None:
            coefficients[column] = (coefficients.get(column, 0) + coeff) % q
    assert total == point(archive["full_sum"])
    assert curve.scalar(total, 4) == target_r
    torsion = curve.add(total, curve.neg(shifted))
    assert torsion in TORSION and curve.add(shifted, torsion) == total
    row = sorted((column, value) for column, value in coefficients.items() if value)
    evaluated = None
    transported = None
    inverse_phase = pow(pow(lam, phase, q), -1, q)
    for column, value in row:
        evaluated = curve.add(evaluated, curve.scalar(column, value))
        transported = curve.add(
            transported, curve.scalar(column, value * inverse_phase % q))
    assert evaluated == target_r
    assert transported == curve.scalar(qpoint, 4)
    return {
        "status": "hit",
        "attempt_position": position,
        "phase": phase,
        "torsion_index": TORSION.index(torsion),
        "witness_indices": indices,
        "full_sum": pjson(total),
        "row": [{"column": pjson(column), "coefficient": value}
                for column, value in row],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    args.out.mkdir(parents=True)
    manifest = json.loads(INPUT.read_text())
    q, lam = manifest["q"], manifest["lambda"]
    assert manifest["field_polynomial"] == 0x80027 and manifest["n"] == 19
    assert manifest["phase_orders"]["spread"] == [(5 * j) % 19 for j in range(19)]
    assert manifest["phase_orders"]["contiguous_control"] == list(range(19))
    assert manifest["caps"] == [1, 4, 8, 19]
    assert sha(ARCHIVE) == "39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c"
    mod = arithmetic()
    field = mod.Field(19, [0, 1, 2, 5])
    field.rabin_prime_degree()
    curve = mod.Curve(field)
    h = tuple(manifest["generator_h"])
    assert mod.source_group_order(19) == 4 * q and mod.is_prime_by_trial(q)
    assert curve.on_curve(h) and curve.scalar(h, q) is None
    assert (lam * lam + lam + 2) % q == 0
    assert lam != 1 and pow(lam, 19, q) == 1
    assert curve.tau(h) == curve.scalar(h, lam)
    assert all(curve.on_curve(t) and curve.scalar(t, 4) is None for t in TORSION)
    factors, rows = histogram(curve, manifest)
    terms, columns = column_terms(curve, factors, q, lam)
    point_index = exact_index(curve, h, q)
    support = bytearray(q)
    inv4 = pow(4, -1, q)
    for projected in rows:
        assert projected in point_index
        support[point_index[projected] * inv4 % q] = 1
    assert sum(support) == manifest["reference_support"] and support[0]
    powers = [pow(lam, s, q) for s in range(19)]
    points = json.loads(POINTS.read_text())
    assert points["q"] == q and points["generator"] == list(h)
    assert len(points["targets"]) == manifest["sample"]["case_count"] == 64
    outcomes = {}
    first_arrays = {}
    for name, order in manifest["phase_orders"].items():
        assert sorted(order) == list(range(19)) and order[0] == 0
        first = first_hits(support, powers, order)
        filename = f"first_hit_{name}.u8"
        (args.out / filename).write_bytes(first)
        first_arrays[name] = first
        certs = []
        for i, case in enumerate(points["targets"]):
            assert case["case_id"] == f"ho-{i:03d}"
            qpoint = tuple(case["point"])
            assert curve.on_curve(qpoint) and qpoint in point_index
            j = point_index[qpoint]
            certs.append({"case_id": case["case_id"],
                          "point": case["point"],
                          "target_index": j,
                          **certificate(curve, qpoint, first[j], order,
                                        factors, rows, terms, q, lam)})
        outcomes[name] = {
            "order": order,
            "first_hit_file": filename,
            "first_hit_sha256": sha(args.out / filename),
            "caps": cap_summary(first, order, manifest["caps"], q),
            "fixed_case_certificates": certs,
        }
    assert outcomes["spread"]["caps"]["19"]["support"] == (
        outcomes["contiguous_control"]["caps"]["19"]["support"])
    spread = outcomes["spread"]["caps"]
    gate = manifest["material_screen"]
    cap4 = spread["4"]
    cap8 = spread["8"]
    pass4 = (cap4["support"] >= gate["primary_cap4_dominance"]["minimum_support"]
             and cap4["oracle_probes"] <= gate["primary_cap4_dominance"]["maximum_oracle_probes"])
    pass8 = (cap8["support"] >= gate["primary_cap8_fallback"]["minimum_support"]
             and cap8["oracle_probes"] * gate["primary_cap8_fallback"]["maximum_oracle_probes_denominator"]
             <= gate["primary_cap8_fallback"]["maximum_oracle_probes_numerator"] * q)
    save(args.out / "outcome.json", {
        "schema": "ecc2k130_rotated_phase_selector_outcome_v1",
        "domain": manifest["domain"],
        "q": q,
        "archive_sha256": sha(ARCHIVE),
        "point_only_sha256": sha(POINTS),
        "factor_sizes": [len(slot) for slot in factors],
        "canonical_nonzero_columns": [pjson(column) for column in columns],
        "single_phase_support": sum(support),
        "orders": outcomes,
        "primary_screen": {"cap4_dominance": pass4, "cap8_fallback": pass8,
                           "priority_pass": pass4 or pass8},
        "field_operations": dict(field.operations),
        "curve_operations": dict(curve.operations),
    })


if __name__ == "__main__":
    main()
