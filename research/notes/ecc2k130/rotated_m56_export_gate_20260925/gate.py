#!/usr/bin/env python3
"""Exact rational-point semantics and capacity gate for frozen rotated PDP targets.

This intentionally does not export or solve a Semaev S6/S7 polynomial. Its
finite-domain reference is an oracle for checking a later algebraic exporter.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import signal
import sys
import tarfile
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS = NOTES / "rotated_pdp_corpus_20260925"
PARENT = NOTES / "rotated_subspace_support_20260925"
ARCHIVE = CORPUS / "evidence/raw.tar.gz"
ARMS = {
    "n13-m5": {"n": 13, "poly": 0x201b, "m": 5, "d": 2, "beta": 3},
    "n19-m6": {"n": 19, "poly": 0x80027, "m": 6, "d": 2, "beta": 3},
}
TORSION = [None, (0, 1), (1, 0), (1, 1)]
WALL_CAP = 600
RSS_CAP = 512 * 1024 * 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canon(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def point(raw):
    return None if raw is None else tuple(raw)


def pjson(p):
    return None if p is None else list(p)


def load_verified_source():
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    inputs = HERE / "INPUTS.json"
    assert sha(inputs) == frozen["inputs_sha256"]
    assert sha(Path(__file__)) == frozen["gate_sha256"]
    assert sha(ARCHIVE) == frozen["archive_sha256"]
    assert sha(CORPUS / "FROZEN.json") == frozen["corpus_frozen_sha256"]
    assert sha(CORPUS / "verify.py") == frozen["corpus_verify_sha256"]
    assert sha(PARENT / "verify.py") == frozen["parent_verify_sha256"]
    manifest = json.loads(inputs.read_text())
    assert manifest["corpus_archive_sha256"] == frozen["archive_sha256"]
    spec = importlib.util.spec_from_file_location("frozen_corpus_verify", CORPUS / "verify.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return manifest, mod


def archive_json(tar: tarfile.TarFile, arm: str, name: str, expected: str | None = None):
    member = tar.getmember(f"raw/{arm}/{name}")
    assert member.isfile() and member.size <= 10_000_000
    contents = tar.extractfile(member).read()
    if expected is not None:
        assert hashlib.sha256(contents).hexdigest() == expected
    return json.loads(contents)


def coordinates(f, beta: int, m: int, d: int):
    conjugates = []
    value = beta
    for _ in range(f.n):
        conjugates.append(value)
        value = f.square(value)
    assert value == beta and f.trace(beta) == 1
    maps = []
    for i in range(m):
        basis = [conjugates[m * j + i] for j in range(d)]
        xs = {}
        for bits in range(1 << d):
            x = 0
            for j, b in enumerate(basis):
                if bits >> j & 1:
                    x ^= b
            assert x not in xs
            xs[x] = bits
            assert f.trace(x) == (bits.bit_count() & 1)
        maps.append(xs)
    assert all(set(maps[i]) & set(maps[j]) == {0}
               for i in range(m) for j in range(i + 1, m))
    return maps


def s3(f, a: int, b: int, c: int) -> int:
    ab = f.mul(a, b)
    return f.mul(f.square(a ^ b), f.square(c)) ^ f.mul(ab, c) ^ f.square(ab) ^ 1


def add_case(p, q) -> str:
    if p is None or q is None:
        return "identity"
    if p[0] != q[0]:
        return "ordinary"
    if p[1] ^ q[1] == p[0]:
        return "inverse_to_infinity"
    assert p == q and p[0] != 0
    return "doubling"


def trace_point(f, p, trace_mask: int) -> int:
    return 0 if p is None else (p[0] & trace_mask).bit_count() & 1


def enumerate_policy(curve, factors, coordinate_maps, target_points, trace_mask,
                     chain: bool):
    f = curve.f
    full = Counter()
    matching = [Counter() for _ in target_points]
    x_masks = [set() for _ in target_points]
    point_targets = {}
    for index, target in enumerate(target_points):
        point_targets.setdefault(target, []).append(index)
    branches = Counter()
    terminal = Counter()
    parity = Counter()
    symmetry_counterexample = None
    total = 0
    for choice in itertools.product(*factors):
        total += 1
        mask_tuple = tuple(coordinate_maps[i][p[0]] for i, p in enumerate(choice))
        bit_parity = sum(v.bit_count() for v in mask_tuple) & 1
        acc = choice[0]
        for p in choice[1:]:
            prior = acc
            case = add_case(prior, p)
            acc = curve.add(prior, p)
            if chain:
                branches[case] += 1
                if prior is not None and acc is not None:
                    assert s3(f, prior[0], p[0], acc[0]) == 0
                elif prior is not None:
                    assert prior[0] == p[0] and acc is None
                else:
                    assert acc == p
        assert bit_parity == trace_point(f, acc, trace_mask)
        parity[bit_parity] += 1
        full[acc] += 1
        for index in point_targets.get(acc, ()):
            matching[index][mask_tuple] += 1
            x_masks[index].add(mask_tuple)
            if chain:
                prior = None
                for p in choice[:-1]:
                    prior = curve.add(prior, p)
                if prior is None:
                    assert choice[-1] == acc
                    terminal["identity_prefix"] += 1
                else:
                    assert s3(f, prior[0], choice[-1][0], acc[0]) == 0
                    terminal["s3_zero"] += 1
        if symmetry_counterexample is None and choice[0][0] and choice[1][0]:
            if choice[0][0] not in coordinate_maps[1] and choice[1][0] not in coordinate_maps[0]:
                symmetry_counterexample = [pjson(choice[0]), pjson(choice[1])]
    projected = {curve.scalar(p, 4) for p in full}
    assert sum(full.values()) == total == math.prod(map(len, factors))
    return {
        "labelled_tuples": total,
        "distinct_full_sums": len(full),
        "distinct_projected_sums": len(projected),
        "parity_counts": dict(sorted(parity.items())),
        "addition_branches": dict(sorted(branches.items())),
        "hit_terminal_branches": dict(sorted(terminal.items())),
        "symmetry_counterexample": symmetry_counterexample,
        "target_full_counts": [sum(c.values()) for c in matching],
        "target_x_mask_sets": [sorted([list(t) for t in masks]) for masks in x_masks],
    }


def run_arm(arm: str):
    start = time.perf_counter()
    cpu = time.process_time()
    def expired(_signum, _frame):
        raise TimeoutError(f"{arm}: {WALL_CAP}s wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, WALL_CAP)
    try:
        manifest, corpus = load_verified_source()
        cfg = ARMS[arm]
        f = corpus.parent_verify.GF(cfg["n"], cfg["poly"])
        curve = corpus.parent_verify.E(f)
        coordinate_maps = coordinates(f, cfg["beta"], cfg["m"], cfg["d"])
        trace_mask = f.trace_mask()
        assert all(f.trace(1 << j) == ((trace_mask >> j) & 1) for j in range(f.n))
        assert [trace_point(f, t, trace_mask) for t in TORSION] == [0, 0, 1, 1]
        assert all(curve.on(t) and curve.scalar(t, 4) is None for t in TORSION)
        with tarfile.open(ARCHIVE, "r:gz") as tar:
            arm_input = manifest["arms"][arm]
            factors_raw = archive_json(tar, arm, "factors.json", arm_input["factor_file_sha256"])
            targets_raw = archive_json(tar, arm, "targets.json", arm_input["target_file_sha256"])
            summary = archive_json(tar, arm, "summary.json")
        assert len(targets_raw) == 8 and len(arm_input["targets"]) == 8
        for row, frozen in zip(targets_raw, arm_input["targets"]):
            assert {k:row[k] for k in frozen} == frozen
        factors = [[point(p) for p in slot] for slot in factors_raw]
        assert factors == corpus.factor_points(curve, cfg["beta"], cfg["m"], cfg["d"])
        assert all({p[0] for p in slot} == set(coordinate_maps[i])
                   for i, slot in enumerate(factors))
        assert all(curve.on(p) for slot in factors for p in slot)
        assert all(len(slot) == (5 if arm == "n13-m5" else 7) for slot in factors)
        target_points = []
        expected_parities = []
        for row in targets_raw:
            q = point(row["Q"])
            assert curve.scalar(q, 4) == point(row["R"])
            assert trace_point(f, q, trace_mask) == 0
            for t in TORSION:
                target = curve.add(q, t)
                assert target is not None
                target_points.append(target)
                expected_parities.append(trace_point(f, t, trace_mask))
                assert trace_point(f, target, trace_mask) == expected_parities[-1]
        rotated = enumerate_policy(curve, factors, coordinate_maps, target_points, trace_mask, True)
        assert rotated["distinct_full_sums"] == summary["distinct_full_sums"]
        assert rotated["distinct_projected_sums"] == summary["distinct_projected_sums"]
        expected_counts = [c for row in targets_raw for c in row["coset_multiplicities"]]
        assert rotated["target_full_counts"] == expected_counts
        assert rotated["symmetry_counterexample"] is not None
        for index, x_masks in enumerate(rotated["target_x_mask_sets"]):
            for mask_tuple in x_masks:
                assert sum(v.bit_count() for v in mask_tuple) & 1 == expected_parities[index]
        repeated_factors = [factors[0]] * cfg["m"]
        repeated_maps = [coordinate_maps[0]] * cfg["m"]
        repeated = enumerate_policy(curve, repeated_factors, repeated_maps,
                                    target_points, trace_mask, False)
        result = {
            "protocol": manifest["domain"], "arm": arm,
            "corpus_commit": "c77767c4a653734f428e110cca29985d721476b2",
            "input_sha256": sha(HERE / "INPUTS.json"),
            "source_sha256": sha(Path(__file__)),
            "archive_sha256": sha(ARCHIVE),
            "trace_mask": trace_mask,
            "target_trace_parities": expected_parities,
            "layout": {
                "direct_summand_bits": cfg["m"] * cfg["d"],
                "chain_intermediate_bits": (cfg["m"] - 2) * cfg["n"],
                "chain_total_bits_before_exception_aux": cfg["m"] * cfg["d"] + (cfg["m"] - 2) * cfg["n"],
                "chain_s3_field_equations": cfg["m"] - 1,
                "chain_boolean_rows_before_exception_aux": (cfg["m"] - 1) * cfg["n"],
                "native_u64_chain_admitted": cfg["m"] * cfg["d"] + (cfg["m"] - 2) * cfg["n"] <= 64,
            },
            "rotated": rotated, "repeated_f0_semantic_control": repeated,
            "wall_seconds": time.perf_counter() - start,
            "cpu_seconds": time.process_time() - cpu,
            "peak_rss_bytes": peak_rss(),
            "field_operations": dict(sorted(f.ops.items())),
            "curve_operations": dict(sorted(curve.ops.items())),
        }
        assert result["wall_seconds"] <= WALL_CAP and result["peak_rss_bytes"] <= RSS_CAP
        return result
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def self_test():
    spec = importlib.util.spec_from_file_location("frozen_corpus_verify", CORPUS / "verify.py")
    assert spec is not None and spec.loader is not None
    corpus = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(corpus)
    f = corpus.parent_verify.GF(13, 0x201b)
    curve = corpus.parent_verify.E(f)
    p = (0, 1)
    assert curve.on(p) and curve.add(p, p) is None
    assert add_case(p, p) == "inverse_to_infinity"
    assert s3(f, 0, 0, 0) == 1  # no affine triple of 2-torsion points sums to O
    assert f.trace(1) == 1
    print("self-test PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=sorted(ARMS))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert args.arm is None and args.out is None
        self_test()
        return
    assert args.arm and args.out
    try:
        result = run_arm(args.arm)
        args.out.write_text(canon(result))
        print(canon({"arm": args.arm, "decision": "PASS", "output_sha256": sha(args.out)}), end="")
    except Exception as error:
        failure = {"arm": args.arm, "decision": "FAIL", "error": repr(error)}
        args.out.write_text(canon(failure))
        print(canon(failure), end="", file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
