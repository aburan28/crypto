#!/usr/bin/env python3
"""Independent bit-polynomial replay of the two frozen n131 slot maps."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sqr_mod(value: int, modulus: int) -> int:
    raw = 0
    while value:
        low = value & -value
        raw ^= 1 << (2 * (low.bit_length() - 1))
        value ^= low
    n = modulus.bit_length() - 1
    while raw.bit_length() > n:
        raw ^= modulus << (raw.bit_length() - n - 1)
    return raw


def rank(values: list[int]) -> int:
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def masks(dimension: int):
    return {"zero": 0, "lowest": 1, "second": 2,
            "lowest_two": 3, "highest": 1 << (dimension - 1),
            "all": (1 << dimension) - 1}


def projected_x(basis: tuple[int, ...], mask: int) -> int:
    x = 0
    for j, value in enumerate(basis):
        if (mask >> j) & 1:
            x ^= value
    return x


def sha_json(value) -> str:
    data = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(data).hexdigest()


def replay() -> dict:
    from capacity import slot_bases_n131

    spec = json.loads((HERE / "INPUT.json").read_text())
    assert spec["field_degree"] == 131 and spec["normal_beta"] == 3
    modulus = int(spec["field_modulus_hex"], 16)
    assert modulus.bit_length() == 132
    conjugates = []
    value = 3
    for _ in range(131):
        conjugates.append(value)
        value = sqr_mod(value, modulus)
    assert value == 3 and rank(conjugates) == 131
    trace = 0
    for value in conjugates:
        trace ^= value
    assert trace == 1
    results = {}
    for arm in spec["arms"]:
        name, dims = arm["name"], arm["dimensions"]
        bases = tuple(tuple(conjugates[10 * j + i] for j in range(d))
                      for i, d in enumerate(dims))
        assert bases == slot_bases_n131(name)
        assert rank([x for basis in bases for x in basis]) == arm["expected_combined_coordinate_rank"]
        assert all(rank(list(basis)) == len(basis) for basis in bases)
        # Check each rotated slot against an independently squared normalized
        # basis, rather than only checking equality to the same index formula.
        for i, basis in enumerate(bases):
            for j, actual in enumerate(basis):
                reference = conjugates[10 * j]
                for _ in range(i):
                    reference = sqr_mod(reference, modulus)
                assert actual == reference
        fixed = []
        for i, basis in enumerate(bases):
            for label, mask in masks(len(basis)).items():
                x = projected_x(basis, mask)
                assert 0 <= x < (1 << 131)
                assert (mask == 0) == (x == 0)
                fixed.append([i, label, mask, hex(x)])
        results[name] = {"slot_dimensions": dims,
                         "combined_coordinate_rank": arm["expected_combined_coordinate_rank"],
                         "basis_sha256": sha_json([[hex(x) for x in slot] for slot in bases]),
                         "fixed_mask_x_sha256": sha_json(fixed),
                         "fixed_mask_count": len(fixed)}
    return {"schema": "ecc2k130-m10-basis-replay-v1",
            "status": "PASS", "normal_basis_rank": 131,
            "normal_beta_trace": 1, "arms": results}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True, separators=(",", ":")))
