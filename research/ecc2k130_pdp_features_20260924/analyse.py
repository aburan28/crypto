#!/usr/bin/env python3
"""Reproduce a toy two-summand PDP feature screen with grouped holdouts."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import sys
import time

SOURCE_SHA = "5259a42b3613e315835546cf3936a4e1dc77c3a7560f3da36dde4e1b233c26d5"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def order_key(prefix: str, value: str) -> bytes:
    return hashlib.sha256((prefix + value).encode("utf-8")).digest()


def orbit_key(field, x: int) -> int:
    assert x != 0
    seen = set()
    while x not in seen:
        seen.add(x)
        x = field.sqr(x)
    assert len(seen) <= field.deg, "orbit length must divide the field degree"
    return min(seen)


class Counts:
    def __init__(self):
        self.mul = self.sqr = self.row_xor = 0

    def as_dict(self):
        return {"field_mul": self.mul, "field_sqr": self.sqr, "row_xor": self.row_xor}


def instrument(field, probe, count: Counts):
    original_mul, original_sqr = field.mul, field.sqr
    original_echelon, original_contains = probe.echelon, probe.contains

    def mul(a, b):
        count.mul += 1
        return original_mul(a, b)

    def sqr(a):
        count.sqr += 1
        return original_sqr(a)

    def echelon(values):
        pivots = {}
        for value in values:
            while value:
                p = value.bit_length() - 1
                if p not in pivots:
                    pivots[p] = value
                    break
                value ^= pivots[p]
                count.row_xor += 1
        return pivots

    def contains(pivots, value):
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                return False
            value ^= pivots[p]
            count.row_xor += 1
        return True

    field.mul, field.sqr = mul, sqr
    probe.echelon, probe.contains = echelon, contains

    def restore():
        field.mul, field.sqr = original_mul, original_sqr
        probe.echelon, probe.contains = original_echelon, original_contains

    return restore


def phase(field, probe, function):
    count = Counts()
    restore = instrument(field, probe, count)
    try:
        before = time.process_time_ns()
        result = function()
        cpu_ns = time.process_time_ns() - before
    finally:
        restore()
    return result, count.as_dict(), cpu_ns


def split_sets(rows, bases, fields):
    base_ids = defaultdict(set)
    targets = defaultdict(set)
    for row in rows:
        n = row["n"]
        base_ids[n].add(row["base_id"])
        targets[n].add(orbit_key(fields[n], row["target"][0]))
    held_bases = {}
    held_orbits = {}
    for n in sorted(base_ids):
        ranked_bases = sorted(base_ids[n], key=lambda b: order_key("pdp-base-holdout-v1|", b))
        assert len(ranked_bases) == 5
        held_bases[n] = ranked_bases[0]
        ranked_orbits = sorted(targets[n], key=lambda x: order_key("pdp-orbit-holdout-v1|", f"{n}|{x}"))
        held_orbits[n] = set(ranked_orbits[:math.ceil(len(ranked_orbits) / 3)])
    return held_bases, held_orbits


def aggregate(rows):
    out = {}
    for split in ("train", "base_holdout", "orbit_holdout", "confirmation"):
        cell = [r for r in rows if r["split"] == split]
        hits = sum(r["group_hit"] for r in cell)
        rejected = sum(r["affine_inconsistent"] for r in cell)
        surviving = [r for r in cell if not r["affine_inconsistent"]]
        exact = all(not r["group_hit"] and r["algebraic_roots"] == 0
                    for r in cell if r["affine_inconsistent"])
        counters = {}
        for unit in ("field_mul", "field_sqr", "row_xor", "cpu_ns"):
            baseline = sum(r["fiber_cost"][unit] for r in cell)
            screened = sum(r["profile_cost"][unit] for r in cell) + sum(
                r["fiber_cost"][unit] for r in surviving)
            counters[unit] = {
                "fiber_only": baseline,
                "screen_plus_remaining_fiber": screened,
                "saved": baseline - screened,
                "screen_to_baseline": (screened / baseline) if baseline else None,
            }
        rank_groups = {}
        for in_span in (False, True):
            group = [r for r in surviving if r["target_in_product_span"] == in_span]
            rank_groups[str(in_span).lower()] = {
                "cases": len(group), "group_hits": sum(r["group_hit"] for r in group),
                "hit_rate": (sum(r["group_hit"] for r in group) / len(group)) if group else None,
            }
        out[split] = {
            "cases": len(cell), "group_hits": hits,
            "certified_affine_refutations": rejected,
            "group_hits_retained": sum(r["group_hit"] for r in surviving),
            "exact_refutation_gate": exact,
            "cost": counters,
            "nonrefuted_by_product_span": rank_groups,
        }
    return out


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-repo", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=root / "producer_results.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "preserve existing evidence"
    source_path = args.source_repo / "scripts/ecc2k130_point_decomposition.py"
    assert sha(source_path) == SOURCE_SHA, "source changed; freeze a new experiment"
    producer = load_module(root / "profile_probe.py", "pdp_profile_probe")
    source = load_module(source_path, "pdp_source")
    artifact = json.loads(args.input.read_text())
    assert artifact["source_sha256"] == SOURCE_SHA
    assert artifact["script_sha256"] == sha(root / "profile_probe.py")
    assert artifact["contract_sha256"] == sha(root / "pdp_profile_contract.md")
    assert artifact["checks"]["cases"] == 716
    natural = [r for r in artifact["records"] if r["kind"] == "natural"]
    assert len(natural) == 640
    bases = {b["base_id"]: b for b in artifact["bases"]}
    assert len(bases) == 20
    fields = {n: source.GF2m(n, next(b["field_polynomial"] for b in bases.values()
                                    if b["base_id"].startswith(f"n{n}-")), tables=True)
              for n in (7, 9, 13, 17)}
    held_bases, held_orbits = split_sets(natural, bases, fields)
    observations = []
    group_checks = 0
    for row in natural:
        n, base_id, target = row["n"], row["base_id"], row["target"]
        F = fields[n]
        basis = bases[base_id]["basis"]
        assert target[0] != 0
        canonical_x = orbit_key(F, target[0])
        bh = base_id == held_bases[n]
        oh = canonical_x in held_orbits[n]
        split = ("confirmation" if bh and oh else
                 "base_holdout" if bh else
                 "orbit_holdout" if oh else "train")
        profile_reps, fiber_reps = [], []
        for _ in range(3):
            (features, affine), pc, pns = phase(F, producer, lambda: producer.profile(F, basis, target[0]))
            (roots, xors), fc, fns = phase(F, producer, lambda: producer.roots_fibers(F, basis, target[0]))
            fc["row_xor"] += xors
            assert features == {k: row[k] for k in features}
            assert len(roots) == row["algebraic_roots"]
            assert features["affine_inconsistent"] is False or not roots
            profile_reps.append((pc, pns))
            fiber_reps.append((fc, fns))
        assert profile_reps[0][0] == profile_reps[1][0] == profile_reps[2][0]
        assert fiber_reps[0][0] == fiber_reps[1][0] == fiber_reps[2][0]
        if split == "confirmation":
            curve = source.Koblitz(F)
            points = {x: curve.points_over(x) for x in producer.span(basis)}
            group_count = len(producer.group_roots(curve, points, tuple(target)))
            assert group_count == row["group_valid_x_pairs"]
            group_checks += 1
        pc = {**profile_reps[0][0], "cpu_ns": int(statistics.median(x[1] for x in profile_reps))}
        fc = {**fiber_reps[0][0], "cpu_ns": int(statistics.median(x[1] for x in fiber_reps))}
        observations.append({
            "base_id": base_id, "n": n, "target_index": row["index"],
            "target": target, "canonical_target_x": canonical_x, "split": split,
            "factor_base_points": bases[base_id]["factor_base_points"],
            "affine_inconsistent": features["affine_inconsistent"],
            "target_in_product_span": features["target_in_product_span"],
            "quadratic_rank": features["quadratic_rank"],
            "product_span_dim": features["product_span_dim"],
            "affine_rank": features["affine_rank"],
            "algebraic_roots": len(roots),
            "group_hit": bool(row["group_valid_x_pairs"]),
            "group_valid_x_pairs": row["group_valid_x_pairs"],
            "profile_cost": pc, "fiber_cost": fc,
            "profile_cpu_ns_reps": [x[1] for x in profile_reps],
            "fiber_cpu_ns_reps": [x[1] for x in fiber_reps],
        })
    summary = aggregate(observations)
    assert all(v["exact_refutation_gate"] for v in summary.values())
    result = {
        "schema": "ecc2k130-pdp-feature-holdout-v1",
        "status": "TOY_DIAGNOSTIC",
        "source_sha256": SOURCE_SHA,
        "producer_sha256": sha(root / "profile_probe.py"),
        "producer_contract_sha256": sha(root / "pdp_profile_contract.md"),
        "analysis_sha256": sha(Path(__file__)),
        "input_sha256": sha(args.input),
        "split": {
            "held_base_by_n": held_bases,
            "held_target_x_orbits_by_n": {n: sorted(v) for n, v in held_orbits.items()},
            "confirmation_group_checks": group_checks,
        },
        "summary": summary,
        "observations": observations,
        "limits": [
            "Toy m=2 systems and unequal useful base point counts",
            "n=17 has subgroup order 239; field degree is not DLP hardness",
            "CPU sums are diagnostic paired observations, not a full-DLP speed claim",
            "Per-query feature calculation recomputes base product span; no caching",
            "The rank feature is a priority hypothesis, not a proof of UNSAT",
        ],
    }
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"status": result["status"], "split": result["split"],
                      "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
