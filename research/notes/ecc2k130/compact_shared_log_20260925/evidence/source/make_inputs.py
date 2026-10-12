#!/usr/bin/env python3
"""Freeze disjoint n37/n41 training scalars and public point-only batch inputs."""
from __future__ import annotations

import argparse
from math import comb
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import types

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
PRIOR = HERE.parent / "compact_base_sweep_20260925"
MATH = HERE.parent / "paired_fullrank_20260925/verify.py"
DOMAIN = b"ECC2K-SPARSE-SHARED-LOG-20260925-v1/"
BASE_COMMIT = "fc27150df3238b6863ed5618c721e7fd8b6ce403"
SOURCE_SHAS = {
    "examples/koblitz_s5_sat_instance.rs": "c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09",
    "examples/koblitz_rho_batch_ks.rs": "fedacb54e441979c8c32860e7b5639e43799234741d49677010164564536d2c8",
}
PRIOR_SPEC_SHA = "700ab214498bdd18cc49530e0233ed2f135628124fdb3eea48cbda9e1ef6666b"
PRIOR_LEDGER_SHA = "d20e66aad75e2f68ed7143ce7df7989666db5d370ef8592f7c6646271cebf1c4"
MATH_SHA = "ed804a6bdace125cc41e514d93b1102763581f4feeaad60ba277cace592230bf"
ARMS = {37: {"q": 230603167, "R": 8, "eta_numerator": 125203,
             "generator_archive": "train-n37-R3-c512"},
        41: {"q": 549756390943, "R": 40, "eta_numerator": 10306,
             "generator_archive": "train-n41-R8-c4096"}}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_once(path: Path, data: bytes, check: bool) -> None:
    if check:
        assert path.is_file(), f"missing frozen input: {path}"
    if path.exists():
        assert path.read_bytes() == data, f"frozen input changed: {path}"
    else:
        path.write_bytes(data)


def reference_math():
    assert sha(MATH.read_bytes()) == MATH_SHA
    module = types.ModuleType("pr737_independent_math")
    module.__file__ = str(MATH)
    exec(compile(MATH.read_bytes(), str(MATH), "exec"), module.__dict__)
    return module


def header_for(n: int, archive_name: str) -> dict:
    path = PRIOR / "evidence/runs" / archive_name / "producer.stdout.jsonl.gz"
    with gzip.open(path, "rb") as stream:
        rows = stream.read().splitlines()
    assert len(rows) == 1
    header = json.loads(rows[0])["compact_orbit_base_header"]
    assert header["n"] == n and header["a"] == 0
    return header


def scalar_stream(n: int, q: int, kind: str, block: int, count: int,
                  seen: set[int]) -> tuple[list[int], list[int]]:
    result, counters = [], []
    counter = 0
    while len(result) < count:
        message = (DOMAIN + str(n).encode() + b"/" + kind.encode() + b"/"
                   + str(block).encode() + b"/" + str(counter).encode())
        value = 1 + int.from_bytes(hashlib.sha256(message).digest(), "big") % (q - 1)
        if value not in seen:
            result.append(value)
            counters.append(counter)
            seen.add(value)
        counter += 1
    return result, counters


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="require every frozen output to exist and match")
    args = parser.parse_args()
    assert sha((PRIOR / "input_spec.json").read_bytes()) == PRIOR_SPEC_SHA
    assert sha((PRIOR / "evidence/SHA256SUMS").read_bytes()) == PRIOR_LEDGER_SHA
    for path, expected in SOURCE_SHAS.items():
        source = subprocess.check_output(["git", "show", f"{BASE_COMMIT}:{path}"], cwd=REPO)
        assert sha(source) == expected, path
    math = reference_math()
    prior_spec = json.loads((PRIOR / "input_spec.json").read_bytes())
    spec = {"schema_version": "1.0", "domain": DOMAIN.decode(), "base_commit": BASE_COMMIT,
            "source_sha256": SOURCE_SHAS, "pr737_math_sha256": MATH_SHA,
            "pr747_input_spec_sha256": PRIOR_SPEC_SHA,
            "pr747_evidence_ledger_sha256": PRIOR_LEDGER_SHA,
            "training_count": 128, "blocks": 3, "max_block_length": 32, "arms": {}}
    for n, arm in ARMS.items():
        q, r, eta = arm["q"], arm["R"], arm["eta_numerator"]
        f = 2*n*r
        assert (2*n*(r-1))**3 * 1_000_000 < 6*q*eta <= f**3 * 1_000_000
        header = header_for(n, arm["generator_archive"])
        assert int(header["subgroup_order"]) == q
        curve = math.Curve(n, 0, header["field_modulus_low_terms"])
        generator = math.pt(header["generator"])
        assert curve.is_on_curve(generator) and curve.scalar(generator, q) is None
        prior_scalars = {int(line) for line in (PRIOR / f"training_scalars_n{n}.txt").read_text().splitlines()}
        assert len(prior_scalars) == 4096
        prior_labels = {int(row["scalar_validator_only"])
                        for holdouts in prior_spec["holdouts"].values() for row in holdouts}
        seen = prior_scalars | prior_labels
        train, train_counters = scalar_stream(n, q, "TRAIN", 0, 128, seen)
        train_data = ("".join(f"{value}\n" for value in train)).encode()
        train_name = f"training_scalars_n{n}.txt"
        write_once(HERE / train_name, train_data, args.check)
        arm_spec = {"n": n, "q": q, "R": r, "F": f,
                    "eta_numerator": eta, "eta_denominator": 1_000_000,
                    "generator": list(generator),
                    "field_modulus_low_terms": header["field_modulus_low_terms"],
                    "generator_archive": arm["generator_archive"],
                    "unordered_multiset_count": comb(f+3, 4),
                    "training_file": train_name, "training_sha256": sha(train_data),
                    "training_counters": train_counters, "training_scalars_validator_only": train,
                    "blocks": []}
        for block in range(3):
            labels, counters = scalar_stream(n, q, "Q", block, 32, seen)
            points = []
            for label in labels:
                point = curve.scalar(generator, label)
                assert point is not None and curve.is_on_curve(point)
                assert curve.scalar(point, q) is None
                points.append(list(point))
            block_spec = {"block": block, "counters": counters,
                          "validator_scalars_only": labels, "points": points, "files": {}}
            for length in (8, 32):
                name = f"points_n{n}_b{block}_L{length}.jsonl"
                data = ("".join(json.dumps(point, separators=(",", ":")) + "\n"
                                for point in points[:length])).encode()
                write_once(HERE / name, data, args.check)
                block_spec["files"][str(length)] = {"name": name, "sha256": sha(data)}
            arm_spec["blocks"].append(block_spec)
        assert len(seen) == 4096 + len({x for x in prior_labels if x not in prior_scalars}) + 128 + 96
        spec["arms"][str(n)] = arm_spec
    data = (json.dumps(spec, indent=2, sort_keys=True) + "\n").encode()
    write_once(HERE / "input_spec.json", data, args.check)
    print(json.dumps({"input_spec_sha256": sha(data),
                      "arms": {n: {"R": spec["arms"][str(n)]["R"],
                                    "train_sha256": spec["arms"][str(n)]["training_sha256"]}
                               for n in ARMS}}, sort_keys=True))


if __name__ == "__main__":
    main()
