#!/usr/bin/env python3
"""Independent group-index replay of the fixed n19 portfolio outcome."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import struct
import tarfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS = NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz"
SWEEP = NOTES / "rotated_beta_sweep_20260925/evidence/raw.tar.gz"
INDEPENDENT = NOTES / "rotated_pdp_corpus_20260925/verify.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(value):
    return None if value is None else tuple(value)


def hist(archive, name: str) -> dict:
    stream = archive.extractfile(name)
    assert stream is not None
    rows = [json.loads(line) for line in stream]
    result = {point(row["point"]): row["count"] for row in rows}
    assert len(result) == len(rows) and all(value > 0 for value in result.values())
    return result


def indexed_patterns(manifest: dict) -> Counter:
    spec = importlib.util.spec_from_file_location("portfolio_independent_curve", INDEPENDENT)
    assert spec is not None and spec.loader is not None
    independent = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(independent)
    field = independent.parent_verify.GF(19, manifest["field_poly"])
    independent.rabin_prime_degree(field)
    curve = independent.parent_verify.E(field)
    q = manifest["q"]
    assert independent.source_order(19) == 4 * q and independent.is_prime(q)
    h = tuple(manifest["reference_h"])
    assert curve.on(h) and curve.scalar(h, q) is None
    four_h = curve.scalar(h, 4)
    assert four_h is not None and curve.scalar(four_h, q) is None

    point_counts, index_counts = [], [None]
    with tarfile.open(CORPUS, "r:gz") as corpus, tarfile.open(SWEEP, "r:gz") as sweep:
        point_counts.append(hist(corpus, "raw/n19-m6/projected_histogram.jsonl"))
        for beta in manifest["beta_order"][1:]:
            prefix = f"raw/beta-{beta}/"
            point_counts.append(hist(sweep, prefix + "projected_histogram.jsonl"))
            stream = sweep.extractfile(prefix + "target_counts.u32le")
            assert stream is not None
            data = stream.read()
            assert len(data) == 4 * q
            index_counts.append(struct.unpack("<" + "I" * q, data))
    assert [len(x) for x in point_counts] == manifest["prior_known_single_support"]
    assert all(sum(row.values()) == manifest["full_tuple_count_per_base"] for row in point_counts)

    patterns = Counter()
    seen = set()
    current = None
    for index in range(q):
        assert current not in seen
        seen.add(current)
        counts = [row.get(current, 0) for row in point_counts]
        assert all(counts[arm] == index_counts[arm][index] for arm in range(1, 5))
        pattern = sum(1 << arm for arm, count in enumerate(counts) if count > 0)
        patterns[pattern] += 1
        current = curve.add(current, four_h)
    assert current is None and len(seen) == q
    assert all(set(row) <= seen for row in point_counts)
    return patterns


def replay(patterns: Counter, manifest: dict) -> dict:
    q, betas = manifest["q"], manifest["beta_order"]
    assert sum(patterns.values()) == q and len(betas) == 5
    unions = {}
    for subset in range(1, 32):
        unions[str(subset)] = q - sum(
            count for pattern, count in patterns.items() if pattern & subset == 0)
    assert [unions[str(1 << i)] for i in range(5)] == manifest["prior_known_single_support"]
    assert [unions[str(1 | (1 << i))] for i in range(1, 5)] == manifest["prior_known_beta3_pair_union"]
    prefixes = []
    previous = 0
    for i, beta in enumerate(betas):
        covered = unions[str((1 << (i + 1)) - 1)]
        prefixes.append({"beta": beta, "base_count": i + 1, "support": covered,
                         "newly_supported": covered - previous, "misses": q - covered})
        previous = covered
    probe_numerator = sum(q - unions[str((1 << i) - 1)] for i in range(1, 5)) + q
    threshold = (manifest["support_followup_threshold_num"] * q +
                 manifest["support_followup_threshold_den"] - 1) // manifest["support_followup_threshold_den"]
    eligible = [subset for subset in range(1, 32) if unions[str(subset)] >= threshold]
    if eligible:
        chosen = min(eligible, key=lambda subset:
                     (subset.bit_count(), tuple(i for i in range(5) if subset & (1 << i))))
        chosen_betas = [beta for i, beta in enumerate(betas) if chosen & (1 << i)]
    else:
        chosen, chosen_betas = None, []
    return {"domain": manifest["domain"], "q": q, "beta_order": betas,
            "pattern_counts": {str(i): patterns[i] for i in range(32)},
            "subset_union_counts": unions, "fixed_order_prefixes": prefixes,
            "sequential_membership_oracle_probe_numerator": probe_numerator,
            "sequential_membership_oracle_probe_denominator": q,
            "support_followup_threshold": threshold,
            "support_followup_pass": bool(eligible),
            "smallest_threshold_subset_mask": chosen,
            "smallest_threshold_subset_betas": chosen_betas}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((HERE / "input_manifest.json").read_text())
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for label, path in (("corpus_archive", CORPUS), ("sweep_archive", SWEEP),
                        ("analyze", HERE / "analyze.py"), ("verify", HERE / "verify.py"),
                        ("independent_curve", INDEPENDENT),
                        ("manifest", HERE / "input_manifest.json"),
                        ("protocol", HERE / "PROTOCOL.md")):
        assert sha(path) == frozen[label + "_sha256"], label
    actual = replay(indexed_patterns(manifest), manifest)
    expected = json.loads(args.expected.read_text())
    assert actual == expected
    args.out.write_text(json.dumps({"status": "PASS", "result_sha256": sha(args.expected),
                                    "all_q_points_checked": manifest["q"],
                                    "four_count_arrays_checked": True},
                                   sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
