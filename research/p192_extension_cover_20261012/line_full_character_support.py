"""Exact character-cell census and known-scalar query replay on 141 classes."""

import argparse
import collections
import hashlib
import json
import math
import sys
import time
from pathlib import Path

from cypari2 import Pari

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_adaptive_private_columns import add, cell, mul, neg, on_curve  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_mod(rows, modulus):
    pivots = {}
    for row in rows:
        vector = dict(collections.Counter(row))
        while vector:
            pivot = min(vector)
            value = vector[pivot] % modulus
            if value == 0:
                del vector[pivot]
                continue
            if pivot not in pivots:
                inverse = pow(value, -1, modulus)
                pivots[pivot] = {k: v * inverse % modulus for k, v in vector.items()
                                 if v % modulus}
                break
            factor = value
            for k, v in pivots[pivot].items():
                vector[k] = (vector.get(k, 0) - factor * v) % modulus
                if vector[k] == 0:
                    del vector[k]
    return len(pivots)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    paths = {
        "base": HERE / "type_i_line_guided_collector.json",
        "density": HERE / "type_i_character_density.json",
        "eligibility": HERE / "type_i_line_fiber_eligibility.json",
        "support": HERE / "type_i_line_three_sum_coverage.json",
        "collector": HERE / "type_i_line_all_cell_collector.json",
    }
    data = {name: json.loads(path.read_text()) for name, path in paths.items()}
    p = int(data["base"]["field_prime"])
    r = int(data["base"]["subgroup_order"])
    base = [tuple(map(int, point)) for point in data["base"]["base_points_fp"]]
    assert len(base) == 141 and all(on_curve(point, p) for point in base)
    coeff = [int(value) for value in data["density"]["sextic_coefficients_decimal"]]
    pari = Pari()
    assert int(pari.isprime(p)) == 1
    pari_p = pari(p)
    halves = [neg(mul(pow(2, -1, r), point, p), p) for point in base]

    def chi(point):
        if point is None:
            return "I"
        value = 0
        for coefficient in reversed(coeff):
            value = (value * point[0] + coefficient) % p
        if value == 0:
            return "0"
        return "+" if int(pari.kronecker(value, pari_p)) == 1 else "-"

    def character(h):
        # The archived line collector orders signs as (Q,H), where Q=2H.
        return chi(add(h, h, p)) + chi(h)

    for i, row in enumerate(data["eligibility"]["rows"]):
        half = neg(halves[i], p)
        assert character(half) == row["cell"] == cell(half, coeff, p)[::-1]

    start_ns = time.monotonic_ns()
    counts = collections.Counter()
    first16 = collections.Counter()
    support = {}
    euler_controls = 0
    for i in range(141):
        for j in range(i, 141):
            pair_half = add(halves[i], halves[j], p)
            pair = add(base[i], base[j], p)
            for k in range(j, 141):
                half = add(pair_half, halves[k], p)
                target = add(pair, base[k], p)
                category = character(half)
                counts[category] += 1
                if k < 16:
                    first16[category] += 1
                if sum(counts.values()) % 503 == 0:
                    assert category == cell(half, coeff, p)[::-1]
                    euler_controls += 1
                assert target not in support
                support[target] = (i, j, k)
        if (i + 1) % 16 == 0:
            print("FULL_CHARACTER_PROGRESS", i + 1, 141, len(support), flush=True)

    assert len(support) == math.comb(143, 3) == 477191
    assert len(support) == data["support"]["distinct_three_sum_targets"]
    assert first16 == collections.Counter({"++": 206, "+-": 218,
                                           "-+": 199, "--": 193})
    assert set(counts) == {"++", "+-", "-+", "--"}

    generator = (int("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16),
                 int("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16))
    arms = {}
    for name in ("control", "balanced"):
        queries = data["collector"][name + "_queries"]
        assert len(queries) == 64
        hits = []
        for query in queries:
            a = int(query["scalar"])
            target = neg(mul(2 * a, generator, p), p)
            half = mul(a, generator, p)
            assert add(half, half, p) == neg(target, p)
            assert character(half) == data["collector"]["trial_rows"][query["trial_index"]]["cell"]
            if target in support:
                hits.append({"scalar": str(a), "indices": support[target]})
        assert len(hits) == len(data["collector"][name + "_hits"])
        arms[name] = {"queries": len(queries), "hits": hits,
                      "rank": rank_mod([hit["indices"] for hit in hits], r)}
        assert arms[name]["rank"] == data["collector"][name + "_rank"]

    receipt = {
        "status": "verified_full_141_character_three_sum_support",
        "field_prime": str(p), "subgroup_order": str(r),
        "base_classes": 141, "three_sum_targets": len(support),
        "support_by_Q_H_cell": dict(sorted(counts.items())),
        "first16_by_Q_H_cell": dict(sorted(first16.items())),
        "plus_plus_only_coverage": counts["++"],
        "four_cell_coverage": len(support),
        "euler_character_controls": euler_controls,
        "known_scalar_arms": arms,
        "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
        "source_sha256": sha(Path(__file__)),
        "input_sha256": {name: sha(path) for name, path in paths.items()},
        "boundary": "Exact support and two 64-query known-scalar laws on the fixed 141-class base; larger-base relation yield and complete one-target costs are separate.",
    }
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print("FULL_CHARACTER_SUPPORT_PASS", json.dumps({
        "support_by_cell": receipt["support_by_Q_H_cell"],
        "control_hits": len(arms["control"]["hits"]),
        "balanced_hits": len(arms["balanced"]["hits"]),
        "euler_controls": euler_controls,
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
