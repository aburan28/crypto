#!/usr/bin/env python3
"""Three-sum coverage after requiring fully split degree-three divisors."""

import collections
import hashlib
import itertools
import json
import sys
from pathlib import Path

from verify_adaptive_private_columns import GX, GY, add, cell, key, mul, neg, on_curve


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def old_index(record_index):
    if record_index < 16:
        return record_index
    if record_index < 80:
        return record_index + 16
    if record_index < 96:
        return record_index - 64
    return record_index


def main():
    if len(sys.argv) != 8:
        raise SystemExit("usage: split_support_collector.py BASE GUIDE DIVISORS DENSITY MIXED OUT.json OUT.gp")
    raw = [Path(path).read_bytes() for path in sys.argv[1:6]]
    base, guide, divisors, density, mixed = list(map(json.loads, raw))
    p, r = int(base["field_prime"]), int(mixed["subgroup_order"])
    assert p == int(divisors["field_prime"]) == int(density["p"])
    assert r == int(density["r"])
    assert divisors["status"] == "exact_fixed_origin_divisor_support_audit"
    generator = (GX, GY)
    assert on_curve(generator, p) and mul(r, generator, p) is None
    unknown = [(int(row["x"]), int(row["y"])) for row in base["points"]]
    assert len(unknown) == 16 and all(on_curve(q, p) for q in unknown)
    scalars = [int(value) for value in guide["accepted_scalars_decimal"]]
    assert len(scalars) == 64
    guided = [neg(mul(2 * value, generator, p), p) for value in scalars]
    old = unknown + [neg(q, p) for q in unknown] + guided + [neg(q, p) for q in guided]
    assert len(old) == len({key(q) for q in old}) == 160
    selected = [(record_index, old_index(record_index))
                for record_index, record in enumerate(divisors["records"])
                if record["linear_roots"] == 3]
    assert len(selected) == 20
    indices = [index for _, index in selected]
    assert len({key(old[index]) for index in indices}) == 20
    unknown_indices = [index for index in indices if index < 32]
    known_indices = [index for index in indices if index >= 32]
    assert len(unknown_indices) == 3 and len(known_indices) == 17
    coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]
    half_negative = (r - 1) // 2
    bearing, anchored = {}, {}
    character_cells = {"bearing": collections.Counter(), "anchored": collections.Counter()}
    repeated_formal = 0
    for triple in itertools.combinations_with_replacement(indices, 3):
        q = add(add(old[triple[0]], old[triple[1]], p), old[triple[2]], p)
        assert q is not None
        row = [0] * 16
        scalar = 0
        for index in triple:
            if index < 16:
                row[index] += 1
            elif index < 32:
                row[index - 16] -= 1
            elif index < 96:
                scalar -= 2 * scalars[index - 32]
            else:
                scalar += 2 * scalars[index - 96]
        scalar %= r
        target = bearing if any(row) else anchored
        point_key = key(q)
        if point_key in target:
            repeated_formal += (target[point_key]["row"] != row or
                                target[point_key]["known_scalar"] != scalar)
        target[point_key] = {"point": q, "row": row, "known_scalar": scalar}
        character_cells["bearing" if any(row) else "anchored"][cell(
            mul(half_negative, q, p), coefficients, p)] += 1
    assert len(bearing) + len(anchored) == 1540
    assert not (set(bearing) & set(anchored))
    assert repeated_formal == 0
    assert (len(bearing), len(anchored)) == (571, 969)
    assert character_cells["bearing"]["++"] == 143
    assert character_cells["anchored"]["++"] == 282
    fresh = mixed["independent_accepted_queries"]
    assert len(fresh) == 64
    fresh_bearing = fresh_anchor = 0
    for entry in fresh:
        scalar = int(entry["scalar"])
        pre = mul(scalar, generator, p)
        assert cell(pre, coefficients, p) == "++"
        qkey = key(neg(add(pre, pre, p), p))
        fresh_bearing += qkey in bearing
        fresh_anchor += qkey in anchored
    steered_keys = sorted(point_key for point_key, record in anchored.items()
                          if cell(mul(half_negative, record["point"], p), coefficients, p) == "++")[:64]
    assert len(steered_keys) == 64
    for point_key in steered_keys:
        record = anchored[point_key]
        assert mul(record["known_scalar"], generator, p) == record["point"]
    lower = int(mixed["accepted_nonzero_scalar_lower_bound"])
    assert 143 * (1 << 182) < lower
    gp_lines = [f"P={p};", f"R={r};", f"B={int('64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1', 16)};",
                f"GX={GX};", f"GY={GY};",
                "C=[" + ",".join(map(str, coefficients)) + "];",
                "POINTS=[" + ",".join(f"[{old[index][0]},{old[index][1]}]" for index in indices) + "];",
                "OLDINDEX=[" + ",".join(map(str, indices)) + "];",
                "KNOWN=[" + ",".join(str((-2 * scalars[index - 32] if index < 96 else
                                              2 * scalars[index - 96]) % r) if index >= 32 else "0"
                                    for index in indices) + "];",
                "FRESH=[" + ",".join(entry["scalar"] for entry in fresh) + "];", ""]
    gp_bytes = "\n".join(gp_lines).encode()
    Path(sys.argv[7]).write_bytes(gp_bytes)
    output = {
        "status": "exact_split_divisor_three_sum_collector",
        "field_prime": str(p), "subgroup_order": str(r),
        "selected_old_indices": indices,
        "selected_unknown_signed_classes": len(unknown_indices),
        "selected_known_signed_classes": len(known_indices),
        "all_unordered_triples": 1540,
        "bearing_targets": len(bearing), "anchored_targets": len(anchored),
        "bearing_anchor_overlap": len(set(bearing) & set(anchored)),
        "repeated_targets_with_distinct_formal_vectors": repeated_formal,
        "character_cells": {name: dict(sorted(counts.items())) for name, counts in character_cells.items()},
        "steered_known_scalar_queries": len(steered_keys),
        "steered_raw_hits": len(steered_keys), "steered_unknown_row_rank": 0,
        "fresh_accepted_queries": len(fresh),
        "fresh_bearing_hits": fresh_bearing, "fresh_anchor_hits": fresh_anchor,
        "ideal_uniform_accepted_hit_bound_numerator": 143,
        "ideal_uniform_accepted_hit_bound_denominator": str(lower),
        "hit_bound_lt_two_pow_minus_182": True,
        "input_sha256": list(map(sha, raw)),
        "source_sha256": sha(Path(__file__).read_bytes()),
        "gp_input_sha256": sha(gp_bytes),
    }
    Path(sys.argv[6]).write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print("SPLIT_SUPPORT_COLLECTOR_PASS signed=20 bearing=571 accepted_bearing=143 "
          f"fresh_hits={fresh_bearing} steered_rank=0")


if __name__ == "__main__":
    main()
