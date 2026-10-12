#!/usr/bin/env python3
"""Losslessly retain target outcomes, skip locations and charged aggregates."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_arm(arm):
    masks = [0] * len(arm["records"])
    for i, k, _ in arm["skipped_residuals"]:
        masks[i] |= 1 << k
    cases = []
    for i, r in enumerate(arm["records"]):
        assert masks[i].bit_count() == r["certified_skips"]
        assert r["hit"] == (r["first_witness"] is not None)
        assert r["third_candidates_tested"] == (
            r["first_witness"][2] + 1 if r["hit"] else arm["base_size"])
        cases.append([
            r["first_witness"], r["third_candidates_tested"], masks[i],
            r["pair_lookups"], r["profile_cache_misses"],
            int(r["independent_before_rank_stop"]), r["rank_after"]])
    return {
        key: value for key, value in arm.items()
        if key not in ("records", "skipped_residuals")
    } | {"cases": cases}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("raw", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), "preserve prior evidence"
    data = json.loads(args.raw.read_text())
    attempts = data["parameters"]["attempts"]
    assert len(data["target_coefficients"]) == attempts
    assert all(len(a["records"]) == attempts for a in data["arms"].values())
    keep = {
        key: value for key, value in data.items()
        if key not in ("arms", "source_target_costs", "codomain_target_costs",
                       "platform", "python")
    }
    keep["schema"] = data["schema"] + "-compact-v1"
    keep["raw_full_sha256"] = digest(args.raw)
    keep["compactor_sha256"] = digest(Path(__file__))
    keep["arms"] = {name: compact_arm(arm)
                    for name, arm in data["arms"].items()}
    args.out.write_text(json.dumps(keep, separators=(",", ":"), sort_keys=True) + "\n")
    print(json.dumps({
        "input": str(args.raw), "output": str(args.out),
        "raw_full_sha256": keep["raw_full_sha256"],
        "raw_bytes": args.raw.stat().st_size,
        "compact_bytes": args.out.stat().st_size,
        "attempts": attempts}, indent=2))


if __name__ == "__main__":
    main()
