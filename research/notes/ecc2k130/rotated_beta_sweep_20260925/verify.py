#!/usr/bin/env python3
"""Independent bit-serial/Fermat replay of beta selection and full tuple census."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import signal
import struct
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORPUS = HERE.parent / "rotated_pdp_corpus_20260925"
_spec = importlib.util.spec_from_file_location("corpus_independent_verify", CORPUS / "verify.py")
assert _spec is not None and _spec.loader is not None
pv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pv)

DOMAIN = "ECC2K130-ROTATED-BETA-SWEEP-20260925-v1"
REFERENCE_ARCHIVE = CORPUS / "evidence" / "raw.tar.gz"
Q = 130873
H = (385982, 301867)
CAP_SECONDS = 600
CAP_RSS = 512 * 1024 * 1024
Point = tuple[int, int] | None


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def point(value) -> Point:
    return None if value is None else tuple(value)


def pjson(value: Point):
    return None if value is None else list(value)


def sorted_key(value: Point):
    return (-1, -1) if value is None else value


def load(path: Path):
    return json.loads(path.read_text())


def jsonl(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def reference():
    with tarfile.open(REFERENCE_ARCHIVE, "r:gz") as archive:
        def read(name):
            stream = archive.extractfile("raw/n19-m6/" + name)
            assert stream is not None
            return stream.read().decode("utf8")
        projected = {point(row["point"]): row["count"]
                     for row in (json.loads(line) for line in read("projected_histogram.jsonl").splitlines())}
        targets = json.loads(read("targets.json"))
    assert len(projected) == 62389 and len(targets) == 8
    assert sum(projected.values()) == 7 ** 6
    return projected, targets


def conjugates(field: pv.parent_verify.GF, beta: int) -> tuple[int, ...]:
    out, value = [], beta
    for _ in range(19):
        out.append(value)
        value = field.square(value)
    assert value == beta
    return tuple(out)


def basis(conj: tuple[int, ...]) -> list[list[int]]:
    result = [[conj[i], conj[6 + i]] for i in range(6)]
    assert all(pv.rank(row) == 2 for row in result)
    assert pv.rank([x for row in result for x in row]) == 12
    return result


def factors(curve: pv.parent_verify.E, beta: int):
    field = curve.f
    conj = conjugates(field, beta)
    assert pv.rank(list(conj)) == 19 and field.trace(beta) == 1
    result = []
    for row in basis(conj):
        xs = [0, row[0], row[1], row[0] ^ row[1]]
        assert len(set(xs)) == 4
        result.append(sorted(p for x in xs for p in pv.lifts(curve, x)))
    assert all(result[i + 1] == sorted(curve.tau(p) for p in result[i])
               for i in range(5))
    return result


def columns(curve: pv.parent_verify.E, points) -> tuple[int, int]:
    projected = {curve.scalar(p, 4) for p in points}
    assert None in projected
    quotient = {None if p is None else min(p, pv.neg(p)) for p in projected}
    return len(projected), len(quotient)


def verify_selection(selection: dict) -> dict:
    assert selection["domain"] == DOMAIN and selection["scan_limit"] == 4096
    assert selection["reference_beta"] == 3
    field = pv.parent_verify.GF(19, 0x80027)
    pv.rabin_prime_degree(field)
    curve = pv.parent_verify.E(field)
    base_conj = conjugates(field, 3)
    baseline_orbit = set(base_conj)
    assert pv.rank(list(base_conj)) == 19 and field.trace(3) == 1
    base = [3, base_conj[6]]
    base_points = sorted(p for x in (0, base[0], base[1], base[0] ^ base[1])
                         for p in pv.lifts(curve, x))
    assert len(base_points) == 7
    baseline_signed, baseline_quotient = columns(curve, base_points)
    assert selection["reference_f0_size"] == len(base_points)
    assert selection["reference_projected_signed"] == baseline_signed
    assert selection["reference_projected_sign_quotient"] == baseline_quotient
    attempts, primary = [], []
    seen_beta, selected_orbits = set(), set(baseline_orbit)
    for counter in range(4096):
        payload = f"{DOMAIN}/beta/{counter}".encode("ascii")
        beta = 1 + int.from_bytes(hashlib.sha256(payload).digest(), "big") % ((1 << 19) - 1)
        row = {"counter": counter, "beta": beta}
        if beta in seen_beta:
            row["reason"] = "duplicate_beta"
        elif beta in selected_orbits:
            row["reason"] = "reference_or_selected_orbit"
        else:
            seen_beta.add(beta)
            conj = conjugates(field, beta)
            row["normal_rank"] = pv.rank(list(conj))
            if row["normal_rank"] != 19:
                row["reason"] = "not_normal"
            else:
                row["trace"] = field.trace(beta)
                rows = basis(conj)
                row["slice_ranks"] = [pv.rank(b) for b in rows]
                row["combined_rank"] = pv.rank([v for b in rows for v in b])
                assert row["trace"] == 1 and row["slice_ranks"] == [2] * 6
                assert row["combined_rank"] == 12
                b0 = rows[0]
                points = sorted(p for x in (0, b0[0], b0[1], b0[0] ^ b0[1])
                                for p in pv.lifts(curve, x))
                row["f0_size"] = len(points)
                if len(points) != 7:
                    row["reason"] = "factor_size"
                else:
                    signed, quotient = columns(curve, points)
                    row["projected_signed"] = signed
                    row["projected_sign_quotient"] = quotient
                    row["reason"] = "equal_column" if signed == baseline_signed else "unequal_column"
                    if signed == baseline_signed:
                        primary.append({"counter": counter, "beta": beta,
                                        "projected_signed": signed,
                                        "projected_sign_quotient": quotient,
                                        "selection": "primary"})
                        selected_orbits.update(conj)
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
            selected_orbits.update(conjugates(field, beta))
    assert selection["attempts"] == attempts
    assert selection["selected"] == selected
    assert selection["fallback_used"] == fallback_used
    assert selection["status"] == ("success" if len(selected) >= 3 else "insufficient_candidates")
    assert selection["wall_seconds"] <= 60 and selection["peak_rss_bytes"] <= 128 * 1024 * 1024
    return {"selection_attempts_replayed": len(attempts), "selected_replayed": len(selected),
            "reference_signed_column": baseline_signed,
            "reference_sign_quotient": baseline_quotient}


def verify_arm(beta: int, archive: Path) -> dict:
    started, cpu_started = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"beta {beta}: {CAP_SECONDS}s verifier cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        selection = load(HERE / "selection.json")
        assert beta in [row["beta"] for row in selection["selected"]]
        reference_projected, frozen_targets = reference()
        field = pv.parent_verify.GF(19, 0x80027)
        pv.rabin_prime_degree(field)
        assert pv.source_order(19) == 4 * Q and pv.is_prime(Q)
        curve = pv.parent_verify.E(field)
        assert curve.on(H) and curve.scalar(H, Q) is None
        assert curve.tau(H) == curve.scalar(H, 41811)
        factor_points = factors(curve, beta)
        assert [len(factor) for factor in factor_points] == [7] * 6
        assert load(archive / "factors.json") == [[pjson(p) for p in row] for row in factor_points]
        column_data = [columns(curve, factor) for factor in factor_points]
        full = Counter()
        for choice in itertools.product(*factor_points):
            acc = None
            for p in choice:
                acc = curve.add(acc, p)
            full[acc] += 1
        assert sum(full.values()) == 7 ** 6
        projected = Counter()
        for s, count in full.items():
            projected[curve.scalar(s, 4)] += count
        assert sum(projected.values()) == 7 ** 6
        rows = jsonl(archive / "full_histogram.jsonl")
        assert [point(row["point"]) for row in rows] == sorted(full, key=sorted_key)
        assert {point(row["point"]): row["count"] for row in rows} == dict(full)
        for row in rows:
            assert pv.tuple_sum(curve, factor_points, row["witness_indices"]) == point(row["point"])
        rows = jsonl(archive / "projected_histogram.jsonl")
        assert [point(row["point"]) for row in rows] == sorted(projected, key=sorted_key)
        assert {point(row["point"]): row["count"] for row in rows} == dict(projected)
        for row in rows:
            s = pv.tuple_sum(curve, factor_points, row["witness_indices"])
            assert s == point(row["full_sum"])
            assert curve.scalar(s, 4) == point(row["point"])
        data = (archive / "target_counts.u32le").read_bytes()
        assert len(data) == 4 * Q
        target_counts = struct.unpack("<" + "I" * Q, data)
        current, four_h = None, curve.scalar(H, 4)
        assert four_h is not None
        expected_counts = []
        for _ in range(Q):
            expected_counts.append(projected.get(current, 0))
            current = curve.add(current, four_h)
        assert current is None and tuple(expected_counts) == target_counts
        assert sum(target_counts) == 7 ** 6
        assert sum(count > 0 for count in target_counts) == len(projected)
        saved_targets = load(archive / "fixed_targets.json")
        assert len(saved_targets) == len(frozen_targets) == 8
        torsion = [None, (0, 1), (1, 0), (1, 1)]
        for i, (actual, frozen) in enumerate(zip(saved_targets, frozen_targets)):
            qpoint, r = point(frozen["Q"]), point(frozen["R"])
            assert curve.scalar(qpoint, 4) == r
            cosets = [curve.add(qpoint, t) for t in torsion]
            counts = [full.get(p, 0) for p in cosets]
            assert actual["index"] == i and actual["reference_class"] == frozen["class"]
            assert actual["Q"] == frozen["Q"] and actual["R"] == frozen["R"]
            assert actual["projected_multiplicity"] == projected.get(r, 0) == sum(counts)
            assert actual["coset_multiplicities"] == counts
            for j, p in enumerate(cosets):
                witness = actual["coset_witness_indices"][j]
                if counts[j]:
                    assert witness is not None and pv.tuple_sum(curve, factor_points, witness) == p
                else:
                    assert witness is None
        summary = load(archive / "summary.json")
        reference_set, current_set = set(reference_projected), set(projected)
        union = len(reference_set | current_set)
        energy = sum(count * count for count in projected.values())
        reference_energy = sum(count * count for count in reference_projected.values())
        assert (energy - 7 ** 6) % 2 == 0
        expected = {"beta": beta, "q": Q, "m": 6, "d": 2,
                    "factor_sizes": [7] * 6,
                    "projected_column_signed": [x[0] for x in column_data],
                    "projected_column_sign_quotient": [x[1] for x in column_data],
                    "labelled_tuples": 7 ** 6,
                    "distinct_full_sums": len(full),
                    "distinct_projected_sums": len(projected),
                    "full_collisions": 7 ** 6 - len(full),
                    "projected_collisions": 7 ** 6 - len(projected),
                    "projected_energy": energy,
                    "projected_colliding_pairs": (energy - 7 ** 6) // 2,
                    "effective_support_floor_num": (7 ** 6) ** 2,
                    "effective_support_floor_den": energy,
                    "reference_projected_energy": reference_energy,
                    "reference_effective_support_floor_num": (7 ** 6) ** 2,
                    "reference_effective_support_floor_den": reference_energy,
                    "exact_misses": Q - len(projected),
                    "counting_minimum_misses": Q - 7 ** 6,
                    "reference_support": len(reference_set),
                    "support_intersection": len(reference_set & current_set),
                    "support_union": union,
                    "candidate_only": len(current_set - reference_set),
                    "reference_only": len(reference_set - current_set),
                    "common_misses": Q - union,
                    "fixed_reference_targets_positive": sum(row["projected_multiplicity"] > 0 for row in saved_targets)}
        for name, value in expected.items():
            assert summary[name] == value, (beta, name, summary[name], value)
        assert summary["total_wall_seconds"] <= 300 and summary["peak_rss_bytes"] <= CAP_RSS
        result = {"beta": beta, "labelled_tuples_replayed": 7 ** 6,
                  "full_histogram_replayed": len(full),
                  "projected_histogram_replayed": len(projected),
                  "all_subgroup_targets_replayed": Q,
                  "fixed_point_targets_replayed": len(saved_targets),
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu_started,
                  "peak_rss_bytes": rss(),
                  "operations": {"field": dict(field.ops), "curve": dict(curve.ops)}}
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        return result
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--beta", type=int)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.selection_only:
        started = time.perf_counter()
        def expired(_signal, _frame):
            raise TimeoutError("selection verifier 60s cap")
        signal.signal(signal.SIGALRM, expired)
        signal.setitimer(signal.ITIMER_REAL, 60)
        try:
            result = verify_selection(load(HERE / "selection.json"))
            assert time.perf_counter() - started <= 60 and rss() <= 128 * 1024 * 1024
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        print(json.dumps(result, sort_keys=True, separators=(",", ":")))
        return
    assert args.beta is not None and args.archive is not None and args.out is not None
    result = verify_arm(args.beta, args.archive)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
