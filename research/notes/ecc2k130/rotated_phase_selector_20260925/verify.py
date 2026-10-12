#!/usr/bin/env python3
"""Independent bit-serial/Fermat replay of both exact all-q phase scans."""
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
INDEPENDENT = NOTES / "rotated_pdp_corpus_20260925/verify.py"
POINTS = NOTES / "rotated_four_base_joint_rank_20260925/inputs/point_only.json"
INPUT = HERE / "INPUT.json"
TORSION = (None, (0, 1), (1, 0), (1, 1))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(value):
    return None if value is None else tuple(value)


def pjson(value):
    return None if value is None else list(value)


def neg(value):
    return None if value is None else (value[0], value[0] ^ value[1])


def independent_arithmetic():
    spec = importlib.util.spec_from_file_location("phase_independent_arithmetic", INDEPENDENT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_archive():
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        member = tar.getmember("raw/n19-m6/factors.json")
        assert member.isfile() and member.size < 100_000
        stream = tar.extractfile(member)
        assert stream is not None
        factors = [[point(p) for p in slot] for slot in json.load(stream)]
        member = tar.getmember("raw/n19-m6/projected_histogram.jsonl")
        assert member.isfile() and member.size < 40_000_000
        stream = tar.extractfile(member)
        assert stream is not None
        rows = [json.loads(line) for line in stream]
    return factors, rows


def first_hit_array(support: bytearray, exponents: list[int], order: list[int]) -> bytes:
    # Different loop shape from the producer: remove targets on first hit.
    q = len(support)
    answers = bytearray([255]) * q
    unresolved = list(range(q))
    for position, phase in enumerate(order):
        multiplier = exponents[phase]
        pending = []
        for j in unresolved:
            if support[(j * multiplier) % q]:
                answers[j] = position
            else:
                pending.append(j)
        unresolved = pending
    assert all(0 <= value < 19 or value == 255 for value in answers)
    return bytes(answers)


def cap_counts(first: bytes, order: list[int], caps: list[int], q: int):
    result = {}
    for cap in caps:
        counts = [first.count(position) for position in range(cap)]
        support = sum(counts)
        remaining = q
        probes = 0
        for count in counts:
            probes += remaining
            remaining -= count
        assert remaining == q - support
        assert probes == sum(min(value + 1, cap) if value != 255 else cap
                             for value in first)
        result[str(cap)] = {
            "support": support,
            "nonzero_support": support - int(first[0] < cap),
            "misses": remaining,
            "first_hit_position_counts": counts,
            "phase_at_position": order[:cap],
            "oracle_probes": probes,
            "probe_mean_denominator": q,
        }
    return result


def make_terms(curve, factors, q: int, lam: int):
    terms = []
    universe = set()
    for i, slot in enumerate(factors):
        local = []
        for source in slot:
            original = source
            for _ in range((19 - i) % 19):
                original = curve.tau(original)
            assert original in factors[0]
            projection = curve.scalar(original, 4)
            if projection is None:
                assert original in TORSION and source in TORSION
                local.append((None, 0))
                continue
            other = neg(projection)
            column = projection if projection < other else other
            sign = 1 if column == projection else -1
            coefficient = sign * pow(lam, i, q) % q
            assert curve.scalar(source, 4) == curve.scalar(column, coefficient)
            universe.add(column)
            local.append((column, coefficient))
        terms.append(local)
    base = {column for column, _ in terms[0] if column is not None}
    assert len(base) == 3 and universe == base
    return terms, sorted(base)


def sample_certificate(curve, target, expected_position: int, order: list[int],
                       histogram: dict, factors, terms, q: int, lam: int):
    if expected_position == 255:
        return {"status": "miss", "attempts": 19}
    phase = order[expected_position]
    transformed = target
    for _ in range(phase):
        transformed = curve.tau(transformed)
    assert transformed == curve.scalar(target, pow(lam, phase, q))
    projection = curve.scalar(transformed, 4)
    saved = histogram[projection]
    indices = saved["witness_indices"]
    assert len(indices) == 6 and all(0 <= index < 7 for index in indices)
    total = None
    coefficients = {}
    for i, index in enumerate(indices):
        source = factors[i][index]
        total = curve.add(total, source)
        column, coefficient = terms[i][index]
        if column is not None:
            coefficients[column] = (coefficients.get(column, 0) + coefficient) % q
    assert total == point(saved["full_sum"])
    assert curve.scalar(total, 4) == projection
    torsion = curve.add(total, neg(transformed))
    assert torsion in TORSION and curve.add(transformed, torsion) == total
    row = sorted((column, coefficient) for column, coefficient in coefficients.items()
                 if coefficient)
    unshifted = None
    shifted = None
    back = pow(lam, (19 - phase) % 19, q)
    for column, coefficient in row:
        shifted = curve.add(shifted, curve.scalar(column, coefficient))
        unshifted = curve.add(unshifted, curve.scalar(column, coefficient * back % q))
    assert shifted == projection and unshifted == curve.scalar(target, 4)
    return {
        "status": "hit",
        "attempt_position": expected_position,
        "phase": phase,
        "torsion_index": TORSION.index(torsion),
        "witness_indices": indices,
        "full_sum": pjson(total),
        "row": [{"column": pjson(column), "coefficient": coefficient}
                for column, coefficient in row],
    }


def replay(expected: dict, evidence: Path) -> dict:
    manifest = json.loads(INPUT.read_text())
    q, lam = manifest["q"], manifest["lambda"]
    assert manifest["phase_orders"]["spread"] == [(5 * j) % 19 for j in range(19)]
    assert manifest["phase_orders"]["contiguous_control"] == list(range(19))
    assert manifest["caps"] == [1, 4, 8, 19]
    assert sha(ARCHIVE) == "39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c"
    assert expected["schema"] == "ecc2k130_rotated_phase_selector_outcome_v1"
    assert expected["domain"] == manifest["domain"] and expected["q"] == q
    assert expected["archive_sha256"] == sha(ARCHIVE)
    assert expected["point_only_sha256"] == sha(POINTS)
    independent = independent_arithmetic()
    field = independent.parent_verify.GF(19, manifest["field_polynomial"])
    independent.rabin_prime_degree(field)
    curve = independent.parent_verify.E(field)
    assert independent.source_order(19) == 4 * q and independent.is_prime(q)
    h = tuple(manifest["generator_h"])
    assert curve.on(h) and curve.scalar(h, q) is None
    assert lam != 1 and pow(lam, 19, q) == 1
    assert (lam * lam + lam + 2) % q == 0 and curve.tau(h) == curve.scalar(h, lam)
    assert all(curve.on(t) and curve.scalar(t, 4) is None for t in TORSION)
    archived_factors, archived_rows = read_archive()
    factors = independent.factor_points(curve, manifest["beta"], manifest["m"],
                                        manifest["d"])
    assert factors == archived_factors
    assert [len(slot) for slot in factors] == expected["factor_sizes"] == [7] * 6
    terms, columns = make_terms(curve, factors, q, lam)
    assert expected["canonical_nonzero_columns"] == [pjson(column) for column in columns]
    histogram = {}
    total = 0
    for row in archived_rows:
        projection = point(row["point"])
        assert projection not in histogram and row["count"] > 0
        histogram[projection] = row
        total += row["count"]
    assert len(histogram) == expected["single_phase_support"] == 62389
    assert total == manifest["tuple_count"] and None in histogram
    # Enumerate Q and [4]Q by independent group addition. This avoids the
    # producer's inverse-four conversion from archived projected indices.
    q_index = {}
    support = bytearray(q)
    current_q = None
    current_r = None
    four_h = curve.scalar(h, 4)
    seen_r = set()
    for j in range(q):
        assert current_q not in q_index and current_r not in seen_r
        q_index[current_q] = j
        seen_r.add(current_r)
        support[j] = int(current_r in histogram)
        current_q = curve.add(current_q, h)
        current_r = curve.add(current_r, four_h)
    assert current_q is None and current_r is None
    assert set(histogram) <= seen_r and sum(support) == 62389 and support[0]
    powers = [pow(lam, i, q) for i in range(19)]
    targets = json.loads(POINTS.read_text())
    assert targets["q"] == q and targets["generator"] == list(h)
    assert len(targets["targets"]) == 64
    all_caps = {}
    first_hashes = {}
    semantic_hits = 0
    for name, order in manifest["phase_orders"].items():
        assert sorted(order) == list(range(19))
        first = first_hit_array(support, powers, order)
        row = expected["orders"][name]
        filename = f"first_hit_{name}.u8"
        assert row["order"] == order and row["first_hit_file"] == filename
        archived_first = (evidence / filename).read_bytes()
        assert len(archived_first) == q and archived_first == first
        assert row["first_hit_sha256"] == sha(evidence / filename)
        first_hashes[name] = row["first_hit_sha256"]
        caps = cap_counts(first, order, manifest["caps"], q)
        assert row["caps"] == caps
        all_caps[name] = caps
        actual_certificates = []
        for i, case in enumerate(targets["targets"]):
            assert case["case_id"] == f"ho-{i:03d}"
            target = tuple(case["point"])
            assert curve.on(target) and target in q_index
            j = q_index[target]
            certificate = sample_certificate(curve, target, first[j], order,
                                             histogram, factors, terms, q, lam)
            if certificate["status"] == "hit":
                semantic_hits += 1
            actual_certificates.append({"case_id": case["case_id"],
                                        "point": case["point"],
                                        "target_index": j,
                                        **certificate})
        assert row["fixed_case_certificates"] == actual_certificates
    assert all_caps["spread"]["19"]["support"] == all_caps["contiguous_control"]["19"]["support"]
    gate = manifest["material_screen"]
    cap4 = all_caps["spread"]["4"]
    cap8 = all_caps["spread"]["8"]
    pass4 = (cap4["support"] >= gate["primary_cap4_dominance"]["minimum_support"]
             and cap4["oracle_probes"] <= gate["primary_cap4_dominance"]["maximum_oracle_probes"])
    pass8 = (cap8["support"] >= gate["primary_cap8_fallback"]["minimum_support"]
             and cap8["oracle_probes"] * gate["primary_cap8_fallback"]["maximum_oracle_probes_denominator"]
             <= gate["primary_cap8_fallback"]["maximum_oracle_probes_numerator"] * q)
    assert expected["primary_screen"] == {
        "cap4_dominance": pass4, "cap8_fallback": pass8,
        "priority_pass": pass4 or pass8,
    }
    return {"status": "PASS", "q_checked": q, "orders": list(manifest["phase_orders"]),
            "first_hit_sha256": first_hashes, "fixed_case_rows_checked": 128,
            "fixed_case_hits_checked": semantic_hits,
            "outcome_sha256": sha(evidence / "outcome.json")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    expected = json.loads((args.evidence / "outcome.json").read_text())
    result = replay(expected, args.evidence)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
