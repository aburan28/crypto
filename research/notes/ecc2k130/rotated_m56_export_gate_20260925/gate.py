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
    "n13-m5": {"n": 13, "poly": 0x201b, "m": 5, "d": 2, "beta": 3, "q": 2003},
    "n19-m6": {"n": 19, "poly": 0x80027, "m": 6, "d": 2, "beta": 3, "q": 130873},
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


def symmetry_controls(curve, factors, coordinate_maps, projected_by_full,
                      ordered_full, excluded_r, q):
    """First exact projected labels lost under two naive slot-order rules."""
    projected = set(projected_by_full.values())
    candidates = {}
    for mode in ("field_integer", "local_mask"):
        ordered_r = {projected_by_full[s] for s in ordered_full[mode]}
        candidates[mode] = {r for r in projected - ordered_r
                            if r is not None and r not in excluded_r}
    first_witness = {mode: {} for mode in candidates}
    for choice in itertools.product(*factors):
        x0, x1 = choice[0][0], choice[1][0]
        if x0 == 0 or x1 == 0:
            continue
        masks = tuple(coordinate_maps[i][p[0]] for i, p in enumerate(choice))
        needed_integer = x0 > x1
        needed_mask = masks[0] > masks[1]
        if not needed_integer and not needed_mask:
            continue
        full_sum = None
        for p in choice:
            full_sum = curve.add(full_sum, p)
        r = projected_by_full[full_sum]
        if needed_integer and r in candidates["field_integer"]:
            first_witness["field_integer"].setdefault(r, (choice, masks, full_sum))
        if needed_mask and r in candidates["local_mask"]:
            first_witness["local_mask"].setdefault(r, (choice, masks, full_sum))
    controls = {}
    for mode in candidates:
        selected = None
        if first_witness[mode]:
            r = min(first_witness[mode])
            choice, masks, full_sum = first_witness[mode][r]
            swapped = (choice[1], choice[0], *choice[2:])
            swapped_sum = None
            for p in swapped:
                swapped_sum = curve.add(swapped_sum, p)
            assert swapped_sum == full_sum
            assert choice[1][0] not in coordinate_maps[0]
            assert choice[0][0] not in coordinate_maps[1]
            subgroup_q = curve.scalar(r, pow(4, -1, q))
            assert curve.scalar(subgroup_q, 4) == r
            torsion = curve.add(full_sum, (subgroup_q[0], subgroup_q[0] ^ subgroup_q[1]))
            assert torsion in TORSION and curve.add(subgroup_q, torsion) == full_sum
            selected = {
                "R": pjson(r), "Q": pjson(subgroup_q), "T": pjson(torsion),
                "full_sum": pjson(full_sum),
                "point_tuple": [pjson(p) for p in choice],
                "x_masks": list(masks),
                "swapped_point_tuple": [pjson(p) for p in swapped],
                "swapped_sum": pjson(swapped_sum),
                "swapped_slot_membership": [False, False],
                "ordered_projected_multiplicity": 0,
            }
        controls[mode] = {
            "projected_labels_lost_by_order": len(candidates[mode]),
            "first_nonzero_two_slot_counterexample": selected,
        }
    return controls


def enumerate_policy(curve, factors, coordinate_maps, target_points, trace_mask,
                     chain: bool, q: int | None = None, excluded_r: set | None = None):
    f = curve.f
    full = Counter()
    matching = [Counter() for _ in target_points]
    x_masks = [set() for _ in target_points]
    mask_witnesses = [{} for _ in target_points]
    point_targets = {}
    for index, target in enumerate(target_points):
        point_targets.setdefault(target, []).append(index)
    branches = Counter()
    terminal = Counter()
    parity = Counter()
    ordered_full = {"field_integer": set(), "local_mask": set()}
    total = 0
    seen_masks = set()
    for choice in itertools.product(*factors):
        total += 1
        mask_tuple = tuple(coordinate_maps[i][p[0]] for i, p in enumerate(choice))
        seen_masks.add(mask_tuple)
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
        if chain:
            if choice[0][0] <= choice[1][0]:
                ordered_full["field_integer"].add(acc)
            if mask_tuple[0] <= mask_tuple[1]:
                ordered_full["local_mask"].add(acc)
        for index in point_targets.get(acc, ()):
            matching[index][mask_tuple] += 1
            x_masks[index].add(mask_tuple)
            mask_witnesses[index].setdefault(mask_tuple, choice)
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
    projected_by_full = {p: curve.scalar(p, 4) for p in full}
    projected = set(projected_by_full.values())
    controls = (symmetry_controls(curve, factors, coordinate_maps,
                                  projected_by_full, ordered_full, excluded_r, q)
                if chain else {})
    assert sum(full.values()) == total == math.prod(map(len, factors))
    assert len(seen_masks) == 4 ** len(factors)
    for index, models in enumerate(mask_witnesses):
        assert set(models) == x_masks[index] == set(matching[index])
        for masks, choice in models.items():
            assert all(p in factors[i] and coordinate_maps[i][p[0]] == masks[i]
                       for i, p in enumerate(choice))
            acc = None
            for p in choice:
                acc = curve.add(acc, p)
            assert acc == target_points[index]
    if chain:
        assert sum(branches.values()) == (len(factors) - 1) * total
        assert sum(terminal.values()) == sum(sum(c.values()) for c in matching)
    else:
        assert not branches and not terminal
    return {
        "labelled_tuples": total,
        "distinct_x_mask_tuples": len(seen_masks),
        "distinct_full_sums": len(full),
        "distinct_projected_sums": len(projected),
        "parity_counts": dict(sorted(parity.items())),
        "addition_branches": dict(sorted(branches.items())),
        "hit_terminal_branches": dict(sorted(terminal.items())),
        "symmetry_negative_controls": controls,
        "target_full_counts": [sum(c.values()) for c in matching],
        "target_x_mask_sets": [sorted([list(t) for t in masks]) for masks in x_masks],
        "target_x_mask_models": [
            [{"masks": list(masks), "rational_lift_multiplicity": count,
              "point_witness": [pjson(p) for p in mask_witnesses[index][masks]]}
             for masks, count in sorted(counter.items())]
            for index, counter in enumerate(matching)
        ],
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
        excluded_r = {point(row["R"]) for row in targets_raw}
        rotated = enumerate_policy(curve, factors, coordinate_maps, target_points,
                                   trace_mask, True, cfg["q"], excluded_r)
        assert rotated["distinct_full_sums"] == summary["distinct_full_sums"]
        assert rotated["distinct_projected_sums"] == summary["distinct_projected_sums"]
        expected_counts = [c for row in targets_raw for c in row["coset_multiplicities"]]
        assert rotated["target_full_counts"] == expected_counts
        for j, row in enumerate(targets_raw):
            if row["class"] == "negative":
                assert all(not masks for masks in rotated["target_x_mask_sets"][4*j:4*j+4])
        for index, x_masks in enumerate(rotated["target_x_mask_sets"]):
            for mask_tuple in x_masks:
                assert sum(v.bit_count() for v in mask_tuple) & 1 == expected_parities[index]
        repeated_factors = [factors[0]] * cfg["m"]
        repeated_maps = [coordinate_maps[0]] * cfg["m"]
        repeated = enumerate_policy(curve, repeated_factors, repeated_maps,
                                    target_points, trace_mask, False)
        if arm == "n13-m5":
            assert repeated["distinct_projected_sums"] == 61  # merged #762 control
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
