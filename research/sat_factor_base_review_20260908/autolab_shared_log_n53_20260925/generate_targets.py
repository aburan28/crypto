#!/usr/bin/env python3
"""Freeze three disjoint public point batches without passing scalar labels to solvers."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

N = 53
MASK = (1 << N) - 1
POLY = (1 << N) | (1 << 6) | (1 << 2) | (1 << 1) | 1
ORDER = 21044858204113
GENERATOR = (198217578752339, 7929897206038174)
TRAIN_DOMAIN = b"ECC2K53-COMPACT-COLD-BATCH-20260924-v1/"
DOMAIN = b"ECC2K53-TRUE-SHARED-LOG-20260925-v1/"


def mul(left: int, right: int) -> int:
    assert 0 <= left <= MASK and 0 <= right <= MASK
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> N:
            left ^= POLY
    return result


def square(value: int) -> int:
    return mul(value, value)


def power(base: int, exponent: int) -> int:
    answer = 1
    while exponent:
        if exponent & 1:
            answer = mul(answer, base)
        base = square(base)
        exponent >>= 1
    return answer


def inverse(value: int) -> int:
    assert value
    answer = power(value, (1 << N) - 2)
    assert mul(value, answer) == 1
    return answer


def on_curve(point: tuple[int, int]) -> bool:
    x, y = point
    return square(y) ^ mul(x, y) == mul(square(x), x) ^ 1


def add(left: tuple[int, int] | None, right: tuple[int, int] | None):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 ^ y2 == x1:
            return None
        assert y1 == y2 and x1 != 0
        slope = x1 ^ mul(y1, inverse(x1))
        x3 = square(slope) ^ slope
        y3 = square(x1) ^ mul(slope ^ 1, x3)
        return x3, y3
    slope = mul(y1 ^ y2, inverse(x1 ^ x2))
    x3 = square(slope) ^ slope ^ x1 ^ x2
    y3 = mul(slope, x1 ^ x3) ^ x3 ^ y1
    return x3, y3


def scalar(point: tuple[int, int], number: int):
    result = None
    while number:
        if number & 1:
            result = add(result, point)
        point = add(point, point)
        number >>= 1
    return result


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def training_scalars() -> set[int]:
    scalars = set()
    counter = 0
    while len(scalars) < 512:
        candidate = 1 + int.from_bytes(
            hashlib.sha256(TRAIN_DOMAIN + str(counter).encode()).digest(), "big"
        ) % (ORDER - 1)
        scalars.add(candidate)
        counter += 1
    return scalars


def generate(out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=False)
    assert on_curve(GENERATOR) and scalar(GENERATOR, ORDER) is None
    seen = training_scalars()
    manifest = {
        "schema_version": "1", "n": N, "a": 0, "order": ORDER,
        "generator": GENERATOR, "domain": DOMAIN.decode(),
        "training_domain": TRAIN_DOMAIN.decode(), "training_count": 512,
        "blocks": [],
    }
    for block in range(3):
        values = []
        counter = 0
        while len(values) < 128:
            candidate = 1 + int.from_bytes(
                hashlib.sha256(DOMAIN + str(block).encode() + b"/" + str(counter).encode()).digest(), "big"
            ) % (ORDER - 1)
            counter += 1
            if candidate in seen:
                continue
            seen.add(candidate)
            values.append(candidate)
        points = [scalar(GENERATOR, value) for value in values]
        assert all(point is not None and on_curve(point) for point in points)
        row = {"block": block, "rho_seed": 531320 + block, "candidate_counters_consumed": counter}
        for count in (32, 128):
            name = f"points_b{block}_L{count}.jsonl"
            payload = b"".join(
                (json.dumps(point, separators=(",", ":")) + "\n").encode()
                for point in points[:count]
            )
            (out / name).write_bytes(payload)
            row[f"L{count}_file"] = name
            row[f"L{count}_sha256"] = sha(payload)
        row["validator_scalars"] = values
        row["validator_scalars_sha256"] = sha(
            b"".join(f"{value}\n".encode() for value in values)
        )
        manifest["blocks"].append(row)
    (out / "TARGET_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(args.out)["blocks"], sort_keys=True))
