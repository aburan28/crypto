#!/usr/bin/env python3
"""Independent Fermat-field and GF(2)-linear S3 candidate replay."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import resource
import signal
import sys
import tarfile
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS_VERIFY = NOTES / "rotated_pdp_corpus_20260925/verify.py"
ARCHIVE = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
GATE = NOTES / "rotated_m56_export_gate_20260925"
ARMS = {"n13-m5": (13, 0x201b, 5, 2), "n19-m6": (19, 0x80027, 6, 2)}
TORSION = [None, (0, 1), (1, 0), (1, 1)]
CAP_SECONDS = 600
CAP_RSS = 512 * 1024 * 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def corpus_module():
    spec = importlib.util.spec_from_file_location("s3_corpus_verify", CORPUS_VERIFY)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def inputs(arm: str):
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["verify_sha256"]
    assert sha(CORPUS_VERIFY) == frozen["corpus_verify_sha256"]
    assert sha(ARCHIVE) == frozen["corpus_archive_sha256"]
    assert sha(GATE / "INPUTS.json") == frozen["gate_inputs_sha256"]
    assert sha(GATE / f"evidence/{arm}.json") == frozen["gate_evidence_sha256"][arm]
    manifest = json.loads((GATE / "INPUTS.json").read_text())
    reference = json.loads((GATE / f"evidence/{arm}.json").read_text())["rotated"]
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        factor_bytes = tar.extractfile(f"raw/{arm}/factors.json").read()
        target_bytes = tar.extractfile(f"raw/{arm}/targets.json").read()
    entry = manifest["arms"][arm]
    assert hashlib.sha256(factor_bytes).hexdigest() == entry["factor_file_sha256"]
    assert hashlib.sha256(target_bytes).hexdigest() == entry["target_file_sha256"]
    factors = [[tuple(p) for p in slot] for slot in json.loads(factor_bytes)]
    targets = json.loads(target_bytes)
    assert len(targets) == 8
    return frozen, factors, targets, reference


def s3(f, a: int, b: int, c: int) -> int:
    A = f.square(a ^ b)
    B = f.mul(a, b)
    return f.mul(A, f.square(c)) ^ f.mul(B, c) ^ f.square(B) ^ 1


def artin_schreier_pivots(f):
    """A right inverse of z -> z^2+z on its trace-zero image."""
    pivots = {}
    for j in range(f.n):
        witness = 1 << j
        image = f.square(witness) ^ witness
        while image:
            high = image.bit_length() - 1
            if high in pivots:
                basis_image, basis_witness = pivots[high]
                image ^= basis_image
                witness ^= basis_witness
            else:
                pivots[high] = (image, witness)
                break
    assert len(pivots) == f.n - 1
    return pivots


def solve_artin_schreier(f, pivots, h: int):
    image, witness = h, 0
    while image:
        high = image.bit_length() - 1
        if high not in pivots:
            return None
        basis_image, basis_witness = pivots[high]
        image ^= basis_image
        witness ^= basis_witness
    assert f.square(witness) ^ witness == h
    return witness


def s3_roots(f, pivots, a: int, b: int, cases: Counter, cache: dict):
    pair = (a, b)
    if pair in cache:
        label, answer = cache[pair]
        cases[label] += 1
        return answer
    A = f.square(a ^ b)
    B = f.mul(a, b)
    C = f.square(B) ^ 1
    if A == 0 and B == 0:
        label, answer = "degenerate_no_root", ()
    elif A == 0:
        label, answer = "linear", (f.mul(C, f.inv(B)),)
    elif B == 0:
        value = f.mul(C, f.inv(A))
        root = value
        for _ in range(f.n - 1):
            root = f.square(root)
        label, answer = "unique_square_root", (root,)
    else:
        h = f.mul(f.mul(A, C), f.inv(f.square(B)))
        z = solve_artin_schreier(f, pivots, h)
        if z is None:
            label, answer = "quadratic_trace_one", ()
        else:
            scale = f.mul(B, f.inv(A))
            label = "quadratic_two_roots"
            answer = tuple(sorted((f.mul(scale, z), f.mul(scale, z ^ 1))))
    assert len(answer) == len(set(answer))
    assert all(s3(f, a, b, c) == 0 for c in answer)
    cache[pair] = (label, answer)
    cases[label] += 1
    return answer


def coordinate_maps(f, m: int, d: int):
    orbit = []
    current = 3
    for _ in range(f.n):
        orbit.append(current)
        current = f.square(current)
    assert current == 3
    arrays = []
    inverse = []
    for i in range(m):
        basis = [orbit[m * bit + i] for bit in range(d)]
        values = [0]
        for value in basis:
            values += [x ^ value for x in values]
        assert len(values) == 1 << d and len(set(values)) == 1 << d
        arrays.append(values)
        inverse.append({x: mask for mask, x in enumerate(values)})
    return arrays, inverse


def replay(arm: str, producer: dict):
    freeze, archived_factors, rows, prior = inputs(arm)
    mod = corpus_module()
    n, poly, m, d = ARMS[arm]
    f = mod.parent_verify.GF(n, poly)
    curve = mod.parent_verify.E(f)
    factors = mod.factor_points(curve, 3, m, d)
    assert factors == archived_factors
    arrays, inverse = coordinate_maps(f, m, d)
    targets = [curve.add(tuple(row["Q"]), torsion) for row in rows for torsion in TORSION]
    assert len(targets) == 32 and all(p is not None for p in targets)
    full_index = {p: j for j, p in enumerate(targets)}
    assert len(full_index) == len(targets)
    by_target_x = defaultdict(list)
    for j, point in enumerate(targets):
        by_target_x[point[0]].append(j)
    trace_mask = f.trace_mask()
    assert all(f.trace(x) == (mask.bit_count() & 1)
               for slot in arrays for mask, x in enumerate(slot))
    parity = [(p[0] & trace_mask).bit_count() & 1 for p in targets]
    assert parity == [v for _ in rows for v in (0, 0, 1, 1)]

    # This is an independent full signed-point census; no producer models enter.
    exact = [defaultdict(lambda: {"all": 0, "affine": 0,
                                   "exceptional": 0, "paths": set()})
             for _ in targets]
    target_branches = [Counter() for _ in targets]
    global_branches = Counter()
    point_count = 0
    for choice in itertools.product(*factors):
        point_count += 1
        masks = tuple(inverse[i][p[0]] for i, p in enumerate(choice))
        acc = choice[0]
        prefix = []
        events = []
        for step, p in enumerate(choice[1:], 1):
            if acc is None:
                case = "identity_prefix"
            elif acc[0] != p[0]:
                case = "ordinary"
            elif acc[1] ^ p[1] == acc[0]:
                case = "inverse_to_O"
            else:
                assert acc == p and acc[0] != 0
                case = "doubling"
            events.append(case)
            global_branches[case] += 1
            acc = curve.add(acc, p)
            if step < m - 1:
                prefix.append(None if acc is None else acc[0])
        j = full_index.get(acc)
        if j is None:
            continue
        assert sum(mask.bit_count() for mask in masks) & 1 == parity[j]
        target_branches[j].update(events)
        item = exact[j][masks]
        item["all"] += 1
        if None in prefix:
            item["exceptional"] += 1
        else:
            item["affine"] += 1
            item["paths"].add(tuple(prefix))
    assert point_count == prior["labelled_tuples"]
    assert sum(global_branches.values()) == point_count * (m - 1)
    assert [sum(x["all"] for x in per_target.values()) for per_target in exact] == prior["target_full_counts"]
    assert [[list(mask) for mask in sorted(per_target)] for per_target in exact] == prior["target_x_mask_sets"]

    pivots = artin_schreier_pivots(f)
    cache = {}
    cases = Counter()
    candidates = [[] for _ in targets]
    mask_parity = Counter()
    path_expansions = 0
    mask_domain = [sorted({inverse[i][point[0]] for point in factor})
                   for i, factor in enumerate(factors)]
    for masks in itertools.product(*mask_domain):
        xs = [arrays[i][mask] for i, mask in enumerate(masks)]
        current_parity = sum(mask.bit_count() for mask in masks) & 1
        mask_parity[current_parity] += 1
        paths = [(x,) for x in s3_roots(f, pivots, xs[0], xs[1], cases, cache)]
        path_expansions += len(paths)
        for x in xs[2:-1]:
            next_paths = []
            for path in paths:
                next_paths.extend(path + (root,) for root in s3_roots(
                    f, pivots, path[-1], x, cases, cache))
            paths = next_paths
            path_expansions += len(paths)
        for path in paths:
            for endpoint in s3_roots(f, pivots, path[-1], xs[-1], cases, cache):
                for j in by_target_x.get(endpoint, ()):
                    candidates[j].append((masks, path, current_parity))

    lift_cache = {}
    comparable_rows = []
    for j, point in enumerate(targets):
        rows_out = []
        candidate_masks = set()
        for masks, path, current_parity in sorted(candidates[j]):
            candidate_masks.add(masks)
            for x in path:
                if x not in lift_cache:
                    lift_cache[x] = bool(mod.lifts(curve, x))
            rows_out.append({"masks": list(masks), "prefix_x": list(path),
                             "trace_pass": current_parity == parity[j],
                             "rational_prefix_lifts": all(lift_cache[x] for x in path),
                             "exact_signed_witness": masks in exact[j] and path in exact[j][masks]["paths"]})
        true_masks = set(exact[j])
        affine_masks = {mask for mask, item in exact[j].items() if item["affine"]}
        exception_only = true_masks - affine_masks
        candidate_pairs = {(tuple(x["masks"]), tuple(x["prefix_x"])) for x in rows_out}
        assert all((mask, path) in candidate_pairs
                   for mask, item in exact[j].items() for path in item["paths"])
        candidate_only = candidate_masks - true_masks
        comparable_rows.append({
            "index": j, "Q_index": j // 4, "T_index": j % 4,
            "target": list(point), "target_x": point[0],
            "target_class": rows[j // 4]["class"], "trace_parity": parity[j],
            "all_mask_count": sum(mask_parity.values()),
            "trace_accepted_mask_count": mask_parity[parity[j]],
            "trace_rejected_mask_count": mask_parity[1 ^ parity[j]],
            "candidate_path_count": len(rows_out),
            "candidate_path_after_trace_count": sum(x["trace_pass"] for x in rows_out),
            "candidate_mask_count": len(candidate_masks),
            "candidate_only_mask_count": len(candidate_only),
            "candidate_nonrational_prefix_path_count": sum(not x["rational_prefix_lifts"] for x in rows_out),
            "candidate_rational_prefix_no_signed_witness_count": sum(
                x["rational_prefix_lifts"] and not x["exact_signed_witness"] for x in rows_out),
            "true_point_tuple_count": sum(x["all"] for x in exact[j].values()),
            "true_affine_point_tuple_count": sum(x["affine"] for x in exact[j].values()),
            "true_exceptional_point_tuple_count": sum(x["exceptional"] for x in exact[j].values()),
            "true_mask_count": len(true_masks),
            "true_affine_mask_count": len(affine_masks),
            "true_exceptional_only_mask_count": len(exception_only),
            "target_point_addition_branches": dict(sorted(target_branches[j].items())),
            "candidate_paths": rows_out,
            "true_exceptional_only_masks": [list(x) for x in sorted(exception_only)],
            "candidate_only_masks": [list(x) for x in sorted(candidate_only)],
        })
    comparisons = {
        "field_degree": n, "field_poly": poly,
        "factor_sizes": list(map(len, factors)),
        "rational_point_tuple_count": point_count,
        "rational_x_mask_tuple_count": sum(mask_parity.values()),
        "mask_parity_counts": dict(sorted(mask_parity.items())),
        "all_point_addition_branches": dict(sorted(global_branches.items())),
        "root_equation_cases": dict(sorted(cases.items())),
        "root_calls": sum(cases.values()), "path_expansions": path_expansions,
        "targets": comparable_rows,
    }
    for key, expected in comparisons.items():
        assert producer[key] == expected, (arm, key, producer.get(key), expected)
    assert producer["protocol"] == freeze["domain"] and producer["arm"] == arm
    return {"arm": arm, "decision": "PASS", "compared_fields": list(comparisons),
            "candidate_path_count": sum(len(x) for x in candidates),
            "exact_signed_point_tuples": sum(row["true_point_tuple_count"] for row in comparable_rows),
            "field_operations": dict(sorted(f.ops.items())),
            "curve_operations": dict(sorted(curve.ops.items()))}


def self_test():
    mod = corpus_module()
    f = mod.parent_verify.GF(3, 0b1011)
    pivots = artin_schreier_pivots(f)
    cases = Counter()
    cache = {}
    for a in range(8):
        for b in range(8):
            expected = tuple(c for c in range(8) if s3(f, a, b, c) == 0)
            assert s3_roots(f, pivots, a, b, cases, cache) == expected
    assert sum(cases.values()) == 64
    print("independent S3 root self-test PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=sorted(ARMS))
    parser.add_argument("--producer", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert args.arm is None and args.producer is None and args.out is None
        self_test()
        return
    assert args.arm and args.producer and args.out
    started, cpu = time.perf_counter(), time.process_time()
    def expired(_signum, _frame):
        raise TimeoutError(f"{args.arm}: {CAP_SECONDS}s verifier wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        producer = json.loads(args.producer.read_text())
        result = replay(args.arm, producer)
        result.update({"wall_seconds": time.perf_counter() - started,
                       "cpu_seconds": time.process_time() - cpu,
                       "peak_rss_bytes": rss(),
                       "producer_sha256": sha(args.producer)})
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        save(args.out, result)
    except Exception as error:
        save(args.out, {"arm": args.arm, "decision": "FAIL",
                        "error": repr(error), "wall_seconds": time.perf_counter() - started,
                        "cpu_seconds": time.process_time() - cpu,
                        "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    main()
