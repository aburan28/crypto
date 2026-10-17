#!/usr/bin/env python3
"""Exact multiset support ceilings for the 648 archived primary cross-stream pairs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
OVERLAP = STUDY / "verification/cross-stream-pilot-overlap-20261015/overlap.json"
AUDIT = STUDY / "verification/cross-stream-pilot-overlap-20261015/audit.json"
COMPARISON = STUDY / "verification/cross-stream-pilot-overlap-20261015/comparison.json"
SAVED = HERE / "screen.json"
INPUT_SHA256 = {
    "overlap.json": "f395be7b22c6d5a05bf1bc6c725ebc927420bacbc2408ac8d9169fe6ed0e4ea2",
    "audit.json": "8ee5622cb6da482993d7cd7b2db3e9c5102367dd71805a6f061167758f2c8235",
    "comparison.json": "5a01a01dd06847c6d330ff82c69ad46c9826800da964fe5ad2d696adf8461dfa",
}
ORDER = 2417851639230796216685689
POINTS_PER_SIGNED_ORBIT = 166
SHAPES = ((5, 1), (5, 2), (6, 1), (6, 2))
BINS = ("under_1_percent", "1_to_10_percent", "10_to_50_percent",
        "50_to_100_percent", "vacuous_at_100_percent")


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiset(points: int, count: int) -> int:
    require(points >= 0 and count >= 0, "negative multiset input")
    return math.comb(points + count - 1, count) if points else int(count == 0)


def support(base_orbits: int, residual_orbits: int, total: int, residuals: int) -> int:
    """Count unordered point multisets with repetition in two disjoint domains."""
    require(0 <= residuals <= total, "invalid summand split")
    return (multiset(base_orbits * POINTS_PER_SIGNED_ORBIT, total - residuals)
            * multiset(residual_orbits * POINTS_PER_SIGNED_ORBIT, residuals))


def bucket(count: int) -> str:
    require(count >= 0, "negative support")
    if count * 100 < ORDER:
        return BINS[0]
    if count * 10 < ORDER:
        return BINS[1]
    if count * 2 < ORDER:
        return BINS[2]
    if count < ORDER:
        return BINS[3]
    return BINS[4]


def build(overlap: dict, audit: dict, comparison: dict) -> dict:
    require(overlap.get("schema") == "n83.primary-cross-stream-orbit-overlap/v1"
            and overlap.get("status") == "PASS_27_ARCHIVES_648_EXACT_PAIRS"
            and overlap.get("archive_count") == 27 and overlap.get("pair_count") == 648
            and overlap.get("selected_best_total_runtime") is None,
            "overlap receipt identity differs")
    require(audit.get("schema") ==
            "n83.primary-cross-stream-orbit-overlap-independent-audit/v1"
            and audit.get("status") == "PASS_27_ARCHIVES_648_PAIRS_MINIMUM_POINT_REPLAY"
            and audit.get("archives_checked") == 27
            and audit.get("pairs_checked") == 648
            and audit.get("receipt_sha256") == INPUT_SHA256["overlap.json"],
            "independent overlap audit differs")
    require(comparison.get("schema") == "n83.cross-stream-overlap-prior-comparison/v1"
            and comparison.get("status") == "PASS_EXACT_MATCH_ALL_27_ARCHIVES_AND_648_PAIRS"
            and comparison.get("pair_matches") == 648
            and comparison.get("archive_matches") == 27
            and comparison.get("current_worker_min_two_residual_pairs") == 558,
            "prior census comparison differs")
    archives = overlap["archives"]
    keys = {(row["policy"], row["seed"], row["columns"]) for row in archives}
    require(len(archives) == len(keys) == 27
            and {key[2] for key in keys} == {64, 256, 600},
            "archive identities differ")
    expected = {(base, pool) for base in keys for pool in keys
                if base[:2] != pool[:2]}
    require(len(expected) == 648, "cross-stream universe differs")
    pairs = overlap["pairs"]
    observed = {(tuple(row["base"]), tuple(row["pool"])) for row in pairs}
    require(len(pairs) == len(observed) == 648 and observed == expected,
            "pair universe differs")

    cases = []
    summary = {f"{total}/{residuals}": {key: 0 for key in BINS}
               for total, residuals in SHAPES}
    eligible = 0
    empty = 0
    singleton = 0
    for row in sorted(pairs, key=lambda item: (item["base"], item["pool"])):
        base = tuple(row["base"])
        pool = tuple(row["pool"])
        shared = row["shared_orbits"]
        residual = row["residual_orbits"]
        require(isinstance(shared, int) and isinstance(residual, int)
                and 0 <= shared <= min(base[2], pool[2])
                and residual == pool[2] - shared
                and row["base_only_orbits"] == base[2] - shared
                and row["union_orbits"] == base[2] + pool[2] - shared
                and row["residual_points"] == POINTS_PER_SIGNED_ORBIT * residual,
                "exact residual partition differs")
        empty += int(residual == 0)
        singleton += int(residual == 1)
        eligible += int(residual >= 2)
        shapes = {}
        for total, residuals in SHAPES:
            count = support(base[2], residual, total, residuals)
            label = bucket(count)
            if residual >= 2:
                summary[f"{total}/{residuals}"][label] += 1
            shapes[f"{total}/{residuals}"] = {
                "unordered_multisets": str(count),
                "uniform_whole_subgroup_hit_ceiling": {
                    "numerator": str(min(count, ORDER)), "denominator": str(ORDER),
                },
                "ceiling_bucket": label,
            }
        cases.append({
            "base": list(base), "pool": list(pool),
            "base_points": base[2] * POINTS_PER_SIGNED_ORBIT,
            "residual_orbits": residual,
            "residual_points": residual * POINTS_PER_SIGNED_ORBIT,
            "current_worker_import_admissible": residual >= 2,
            "shapes": shapes,
        })
    require((eligible, empty, singleton) == (558, 33, 57),
            "pilot residual admission counts differ")
    require(all(sum(bins.values()) == eligible for bins in summary.values()),
            "support bucket census incomplete")
    return {
        "schema": "n83.primary-cross-stream-support-ceilings/v1",
        "status": "PASS_EXACT_648_PAIR_SUPPORT_SCREEN",
        "source_sha256": sha256(Path(__file__)),
        "input_sha256": INPUT_SHA256,
        "subgroup_order": str(ORDER),
        "sampling_premise": "uniform whole prime-order subgroup query target",
        "counting_model": "unordered point multisets with repetition in exact disjoint small and residual domains",
        "scope_limit": "counting upper bounds only; no solver yield, graph closure, rank, LA, target log or total runtime measured",
        "pair_count": len(cases), "current_worker_admissible_pairs": eligible,
        "empty_residual_pairs": empty, "singleton_residual_pairs": singleton,
        "admissible_bucket_counts": summary,
        "cases": cases,
        "selected_best_total_runtime": None,
    }


def generate() -> dict:
    paths = {"overlap.json": OVERLAP, "audit.json": AUDIT,
             "comparison.json": COMPARISON}
    for name, path in paths.items():
        require(sha256(path) == INPUT_SHA256[name], f"{name} content changed")
    return build(*(json.loads(paths[name].read_text()) for name in
                   ("overlap.json", "audit.json", "comparison.json")))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = json.dumps(generate(), indent=2, sort_keys=True) + "\n"
    if args.write:
        with SAVED.open("x") as stream:
            stream.write(output)
    else:
        require(SAVED.read_text() == output, "saved support screen differs")
    print(json.dumps({"status":"PASS", "sha256":hashlib.sha256(output.encode()).hexdigest(),
                      "bytes":len(output)}, sort_keys=True))


if __name__ == "__main__":
    main()
