#!/usr/bin/env python3
"""Frozen exact affine recursive-S3 candidate census; see PROTOCOL.md."""
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
ARCHIVE = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
PARENT = NOTES / "rotated_subspace_support_20260925/gate.py"
GATE = NOTES / "rotated_m56_export_gate_20260925"
TORSION = [None, (0, 1), (1, 0), (1, 1)]
ARMS = {
    "n13-m5": (13, 0x201b, [0, 1, 3, 4], 5, 2, 2003),
    "n19-m6": (19, 0x80027, [0, 1, 2, 5], 6, 2, 130873),
}
CAP_SECONDS = 600
CAP_RSS = 512 * 1024 * 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def parent_module():
    spec = importlib.util.spec_from_file_location("s3_parent_producer", PARENT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def frozen_inputs(arm: str):
    freeze = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == freeze["diagnostic_sha256"]
    assert sha(PARENT) == freeze["parent_producer_sha256"]
    assert sha(ARCHIVE) == freeze["corpus_archive_sha256"]
    assert sha(GATE / "INPUTS.json") == freeze["gate_inputs_sha256"]
    assert sha(GATE / "FROZEN.json") == freeze["gate_frozen_sha256"]
    assert sha(GATE / f"evidence/{arm}.json") == freeze["gate_evidence_sha256"][arm]
    manifest = json.loads((GATE / "INPUTS.json").read_text())
    reference = json.loads((GATE / f"evidence/{arm}.json").read_text())["rotated"]
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        factors_content = tar.extractfile(f"raw/{arm}/factors.json").read()
        targets_content = tar.extractfile(f"raw/{arm}/targets.json").read()
    entry = manifest["arms"][arm]
    assert hashlib.sha256(factors_content).hexdigest() == entry["factor_file_sha256"]
    assert hashlib.sha256(targets_content).hexdigest() == entry["target_file_sha256"]
    factors = [[tuple(p) for p in slot] for slot in json.loads(factors_content)]
    targets = json.loads(targets_content)
    assert len(targets) == len(entry["targets"]) == 8
    for row, frozen in zip(targets, entry["targets"]):
        assert {key: row[key] for key in frozen} == frozen
    return freeze, factors, targets, reference


def square_root(f, value: int) -> int:
    root = value
    for _ in range(f.n - 1):
        root = f.square(root)
    assert f.square(root) == value
    return root


def half_trace(f, value: int) -> int:
    assert f.n & 1 and f.trace(value) == 0
    result, term = 0, value
    for _ in range((f.n + 1) // 2):
        result ^= term
        term = f.square(f.square(term))
    assert f.square(result) ^ result == value
    return result


def s3(f, a: int, b: int, c: int) -> int:
    product = f.mul(a, b)
    return f.mul(f.square(a ^ b), f.square(c)) ^ f.mul(product, c) ^ f.square(product) ^ 1


def roots(f, a: int, b: int, cases: Counter) -> tuple[int, ...]:
    product = f.mul(a, b)
    A, B, C = f.square(a ^ b), product, f.square(product) ^ 1
    if A == 0:
        if B == 0:
            cases["degenerate_no_root"] += 1
            out = ()
        else:
            cases["linear"] += 1
            out = (f.mul(C, f.inverse(B)),)
    elif B == 0:
        cases["unique_square_root"] += 1
        out = (square_root(f, f.mul(C, f.inverse(A))),)
    else:
        h = f.mul(f.mul(A, C), f.inverse(f.square(B)))
        if f.trace(h):
            cases["quadratic_trace_one"] += 1
            out = ()
        else:
            cases["quadratic_two_roots"] += 1
            u = half_trace(f, h)
            scale = f.mul(B, f.inverse(A))
            out = tuple(sorted((f.mul(scale, u), f.mul(scale, u ^ 1))))
    assert len(out) == len(set(out))
    assert all(s3(f, a, b, c) == 0 for c in out)
    return out


def lifts(f, x: int) -> tuple[tuple[int, int], ...]:
    if x == 0:
        return ((0, 1),)
    rhs = x ^ f.square(f.inverse(x))
    if f.trace(rhs):
        return ()
    z = half_trace(f, rhs)
    return tuple(sorted(((x, f.mul(x, z)), (x, f.mul(x, z ^ 1)))))


def branch(prior, factor) -> str:
    if prior is None:
        return "identity_prefix"
    if prior[0] != factor[0]:
        return "ordinary"
    if prior[1] ^ factor[1] == prior[0]:
        return "inverse_to_O"
    assert prior == factor and prior[0] != 0
    return "doubling"


def coordinates(f, mod, m: int, d: int):
    conjugates = mod.normal_conjugates(f, 3)
    bases = mod.subspace_basis(conjugates, m, d)
    by_mask = []
    by_x = []
    for basis in bases:
        values = []
        for mask in range(1 << d):
            x = 0
            for bit, element in enumerate(basis):
                if mask >> bit & 1:
                    x ^= element
            values.append(x)
        assert len(set(values)) == 1 << d
        by_mask.append(values)
        by_x.append({x: mask for mask, x in enumerate(values)})
    return by_mask, by_x


def complete_reference(curve, factors, by_x, targets, trace_mask, expected):
    m = len(factors)
    full_index = {p: i for i, p in enumerate(targets)}
    assert len(full_index) == len(targets)
    models = [defaultdict(lambda: {"all": 0, "affine": 0, "exceptional": 0,
                                   "paths": set()}) for _ in targets]
    branches = [Counter() for _ in targets]
    all_branches = Counter()
    point_tuples = 0
    for choice in itertools.product(*factors):
        point_tuples += 1
        masks = tuple(by_x[i][p[0]] for i, p in enumerate(choice))
        acc = choice[0]
        prefixes = []
        events = []
        for j, p in enumerate(choice[1:], start=1):
            case = branch(acc, p)
            events.append(case)
            all_branches[case] += 1
            acc = curve.add(acc, p)
            if j < m - 1:
                prefixes.append(None if acc is None else acc[0])
        index = full_index.get(acc)
        if index is None:
            continue
        assert (sum(mask.bit_count() for mask in masks) & 1) == ((acc[0] & trace_mask).bit_count() & 1)
        branches[index].update(events)
        model = models[index][masks]
        model["all"] += 1
        if any(x is None for x in prefixes):
            model["exceptional"] += 1
        else:
            model["affine"] += 1
            model["paths"].add(tuple(prefixes))
    assert [sum(v["all"] for v in per_target.values()) for per_target in models] == expected["target_full_counts"]
    assert [[list(x) for x in sorted(per_target)] for per_target in models] == expected["target_x_mask_sets"]
    assert sum(all_branches.values()) == (m - 1) * point_tuples
    return models, branches, all_branches, point_tuples


def run_arm(arm: str, out: Path):
    started, cpu_started = time.perf_counter(), time.process_time()
    def expired(_signum, _frame):
        raise TimeoutError(f"{arm}: {CAP_SECONDS}s wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    out.mkdir(parents=True, exist_ok=False)
    try:
        frozen, factors, target_rows, prior_result = frozen_inputs(arm)
        n, poly, low, m, d, _q = ARMS[arm]
        mod = parent_module()
        f = mod.Field(n, low)
        curve = mod.Curve(f)
        assert f.poly == poly
        by_mask, by_x = coordinates(f, mod, m, d)
        assert factors == [sorted(p for x in set(slot) for p in lifts(f, x))
                           for slot in [[p[0] for p in factor] for factor in factors]]
        assert all(curve.on_curve(p) for factor in factors for p in factor)
        target_points = [curve.add(tuple(row["Q"]), torsion) for row in target_rows
                         for torsion in TORSION]
        assert all(p is not None for p in target_points)
        target_x = defaultdict(list)
        for index, p in enumerate(target_points):
            target_x[p[0]].append(index)
        trace_mask = sum(f.trace(1 << bit) << bit for bit in range(n))
        assert all((f.trace(x) == (mask.bit_count() & 1))
                   for slot in by_mask for mask, x in enumerate(slot))
        expected_parity = [(p[0] & trace_mask).bit_count() & 1 for p in target_points]
        assert expected_parity == [v for _ in target_rows for v in (0, 0, 1, 1)]
        models, branches, all_branches, tuple_count = complete_reference(
            curve, factors, by_x, target_points, trace_mask, prior_result)
        assert tuple_count == prior_result["labelled_tuples"]

        mask_domains = [sorted({by_x[i][p[0]] for p in factor})
                        for i, factor in enumerate(factors)]
        cases = Counter()
        candidates = [[] for _ in target_points]
        mask_parity = Counter()
        path_expansions = 0
        for ordinal, masks in enumerate(itertools.product(*mask_domains)):
            xs = [by_mask[i][mask] for i, mask in enumerate(masks)]
            parity = sum(v.bit_count() for v in masks) & 1
            mask_parity[parity] += 1
            paths = [(root,) for root in roots(f, xs[0], xs[1], cases)]
            path_expansions += len(paths)
            for x in xs[2:-1]:
                paths = [path + (root,) for path in paths
                         for root in roots(f, path[-1], x, cases)]
                path_expansions += len(paths)
            for path in paths:
                for terminal_x in roots(f, path[-1], xs[-1], cases):
                    for index in target_x.get(terminal_x, ()):
                        candidates[index].append((masks, path, parity))
            if ordinal % 64 == 0:
                save(out / "progress.json", {"arm": arm, "masks_completed": ordinal + 1,
                                             "path_expansions": path_expansions,
                                             "wall_seconds": time.perf_counter() - started})
        mask_count = sum(mask_parity.values())
        assert mask_count == prior_result["distinct_x_mask_tuples"]
        lift_cache = {}
        rows = []
        for index, point in enumerate(target_points):
            exact = models[index]
            candidate_rows = []
            candidate_masks = set()
            for masks, path, parity in sorted(candidates[index]):
                candidate_masks.add(masks)
                lifts_ok = []
                for x in path:
                    if x not in lift_cache:
                        lift_cache[x] = bool(lifts(f, x))
                    lifts_ok.append(lift_cache[x])
                signed = masks in exact and path in exact[masks]["paths"]
                candidate_rows.append({"masks": list(masks), "prefix_x": list(path),
                                       "trace_pass": parity == expected_parity[index],
                                       "rational_prefix_lifts": all(lifts_ok),
                                       "exact_signed_witness": signed})
            exact_masks = set(exact)
            affine_masks = {mask for mask, data in exact.items() if data["affine"]}
            exceptional_masks = exact_masks - affine_masks
            assert affine_masks <= candidate_masks
            candidate_pairs = {(tuple(row["masks"]), tuple(row["prefix_x"]))
                               for row in candidate_rows}
            assert all((mask, path) in candidate_pairs
                       for mask, data in exact.items() for path in data["paths"])
            candidate_only = candidate_masks - exact_masks
            rows.append({
                "index": index, "Q_index": index // 4, "T_index": index % 4,
                "target": list(point), "target_x": point[0],
                "target_class": target_rows[index // 4]["class"],
                "trace_parity": expected_parity[index],
                "all_mask_count": mask_count,
                "trace_accepted_mask_count": mask_parity[expected_parity[index]],
                "trace_rejected_mask_count": mask_parity[1 ^ expected_parity[index]],
                "candidate_path_count": len(candidate_rows),
                "candidate_path_after_trace_count": sum(row["trace_pass"] for row in candidate_rows),
                "candidate_mask_count": len(candidate_masks),
                "candidate_only_mask_count": len(candidate_only),
                "candidate_nonrational_prefix_path_count": sum(not row["rational_prefix_lifts"] for row in candidate_rows),
                "candidate_rational_prefix_no_signed_witness_count": sum(
                    row["rational_prefix_lifts"] and not row["exact_signed_witness"]
                    for row in candidate_rows),
                "true_point_tuple_count": sum(data["all"] for data in exact.values()),
                "true_affine_point_tuple_count": sum(data["affine"] for data in exact.values()),
                "true_exceptional_point_tuple_count": sum(data["exceptional"] for data in exact.values()),
                "true_mask_count": len(exact_masks),
                "true_affine_mask_count": len(affine_masks),
                "true_exceptional_only_mask_count": len(exceptional_masks),
                "target_point_addition_branches": dict(sorted(branches[index].items())),
                "candidate_paths": candidate_rows,
                "true_exceptional_only_masks": [list(mask) for mask in sorted(exceptional_masks)],
                "candidate_only_masks": [list(mask) for mask in sorted(candidate_only)],
            })
            assert rows[-1]["true_point_tuple_count"] == prior_result["target_full_counts"][index]
        result = {
            "protocol": frozen["domain"], "arm": arm, "field_degree": n,
            "field_poly": poly, "factor_sizes": list(map(len, factors)),
            "rational_point_tuple_count": tuple_count,
            "rational_x_mask_tuple_count": mask_count,
            "mask_parity_counts": dict(sorted(mask_parity.items())),
            "all_point_addition_branches": dict(sorted(all_branches.items())),
            "root_equation_cases": dict(sorted(cases.items())),
            "root_calls": sum(cases.values()), "path_expansions": path_expansions,
            "targets": rows,
            "field_operations": dict(sorted(f.operations.items())),
            "curve_operations": dict(sorted(curve.operations.items())),
            "wall_seconds": time.perf_counter() - started,
            "cpu_seconds": time.process_time() - cpu_started,
            "peak_rss_bytes": rss(),
        }
        assert result["wall_seconds"] <= CAP_SECONDS and result["peak_rss_bytes"] <= CAP_RSS
        save(out / "result.json", result)
        return result
    except Exception as error:
        save(out / "failure.json", {"arm": arm, "error": repr(error),
                                    "wall_seconds": time.perf_counter() - started,
                                    "cpu_seconds": time.process_time() - cpu_started,
                                    "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def self_test():
    mod = parent_module()
    f = mod.Field(3, [0, 1])
    cases = Counter()
    for a in range(8):
        for b in range(8):
            expected = tuple(c for c in range(8) if s3(f, a, b, c) == 0)
            assert roots(f, a, b, cases) == expected
    assert sum(cases.values()) == 64
    print("S3 root self-test PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=sorted(ARMS))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert args.arm is None and args.out is None
        self_test()
    else:
        assert args.arm is not None and args.out is not None
        run_arm(args.arm, args.out)


if __name__ == "__main__":
    main()
