#!/usr/bin/env python3
"""Exact set analysis of the five frozen n19 projected support archives."""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
SWEEP = NOTES / "rotated_beta_sweep_20260925/evidence/raw.tar.gz"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(value):
    return None if value is None else tuple(value)


def archive_sets(manifest: dict) -> list[set]:
    result = []
    with tarfile.open(CORPUS, "r:gz") as corpus, tarfile.open(SWEEP, "r:gz") as sweep:
        for beta in manifest["beta_order"]:
            archive, name = (corpus, "raw/n19-m6/projected_histogram.jsonl") if beta == 3 else (
                sweep, f"raw/beta-{beta}/projected_histogram.jsonl")
            stream = archive.extractfile(name)
            assert stream is not None
            rows = [json.loads(line) for line in stream]
            points = [point(row["point"]) for row in rows]
            counts = [row["count"] for row in rows]
            assert len(points) == len(set(points)) and all(c > 0 for c in counts)
            assert sum(counts) == manifest["full_tuple_count_per_base"]
            result.append(set(points))
    assert [len(x) for x in result] == manifest["prior_known_single_support"]
    assert [len(result[0] | x) for x in result[1:]] == manifest["prior_known_beta3_pair_union"]
    return result


def outcome(sets: list[set], manifest: dict) -> dict:
    betas, q = manifest["beta_order"], manifest["q"]
    assert len(sets) == len(betas) == 5
    patterns = Counter()
    for point_value in set().union(*sets):
        pattern = sum((1 << i) for i, arm in enumerate(sets) if point_value in arm)
        assert 0 < pattern < 32
        patterns[pattern] += 1
    patterns[0] = q - sum(patterns.values())
    assert patterns[0] >= 0 and sum(patterns.values()) == q
    unions = {str(subset): sum(count for pattern, count in patterns.items() if pattern & subset)
              for subset in range(1, 32)}
    assert all(unions[str(1 << i)] == len(sets[i]) for i in range(5))
    prefixes, prior = [], 0
    for i, beta in enumerate(betas):
        mask = (1 << (i + 1)) - 1
        covered = unions[str(mask)]
        prefixes.append({"beta": beta, "base_count": i + 1, "support": covered,
                         "newly_supported": covered - prior, "misses": q - covered})
        prior = covered
    probe_numerator = q + sum(row["misses"] for row in prefixes[:-1])
    threshold = (manifest["support_followup_threshold_num"] * q +
                 manifest["support_followup_threshold_den"] - 1) // manifest["support_followup_threshold_den"]
    eligible = [subset for subset in range(1, 32) if unions[str(subset)] >= threshold]
    chosen = min(eligible, key=lambda x: (x.bit_count(), tuple(i for i in range(5) if x >> i & 1))) if eligible else None
    return {"domain": manifest["domain"], "q": q, "beta_order": betas,
            "pattern_counts": {str(i): patterns[i] for i in range(32)},
            "subset_union_counts": unions, "fixed_order_prefixes": prefixes,
            "sequential_membership_oracle_probe_numerator": probe_numerator,
            "sequential_membership_oracle_probe_denominator": q,
            "support_followup_threshold": threshold,
            "support_followup_pass": bool(eligible),
            "smallest_threshold_subset_mask": chosen,
            "smallest_threshold_subset_betas": [betas[i] for i in range(5) if chosen is not None and chosen >> i & 1]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for label, path in (("corpus_archive", CORPUS), ("sweep_archive", SWEEP),
                        ("analyze", HERE / "analyze.py"), ("verify", HERE / "verify.py"),
                        ("manifest", HERE / "input_manifest.json"),
                        ("protocol", HERE / "PROTOCOL.md")):
        assert sha(path) == frozen[label + "_sha256"], label
    result = outcome(archive_sets(manifest), manifest)
    args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
