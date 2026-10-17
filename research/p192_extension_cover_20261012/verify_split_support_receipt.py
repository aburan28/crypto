#!/usr/bin/env python3
"""Bind the split-divisor collector receipt to its independent PARI replay."""

import hashlib
import json
import re
import sys
from pathlib import Path


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def old_index(index):
    if index < 16:
        return index
    if index < 80:
        return index + 16
    if index < 96:
        return index - 64
    return index


def main():
    if len(sys.argv) != 10:
        raise SystemExit("usage: verify_split_support_receipt.py BASE GUIDE DIVISOR DENSITY MIXED RECEIPT GP_INPUT GP_LOG PRODUCER")
    raw = [Path(name).read_bytes() for name in sys.argv[1:6]]
    divisors = json.loads(raw[2])
    receipt = json.loads(Path(sys.argv[6]).read_bytes())
    gp_input = Path(sys.argv[7]).read_bytes()
    gp_log = Path(sys.argv[8]).read_text()
    assert receipt["status"] == "exact_split_divisor_three_sum_collector"
    assert receipt["input_sha256"] == list(map(sha, raw))
    assert receipt["gp_input_sha256"] == sha(gp_input)
    assert receipt["source_sha256"] == sha(Path(sys.argv[9]).read_bytes())
    selected = [old_index(index) for index, row in enumerate(divisors["records"])
                if row["linear_roots"] == 3]
    assert receipt["selected_old_indices"] == selected
    assert (receipt["selected_unknown_signed_classes"],
            receipt["selected_known_signed_classes"]) == (3, 17)
    assert (receipt["all_unordered_triples"], receipt["bearing_targets"],
            receipt["anchored_targets"], receipt["bearing_anchor_overlap"],
            receipt["repeated_targets_with_distinct_formal_vectors"]) == (1540, 571, 969, 0, 0)
    assert (receipt["fresh_accepted_queries"], receipt["fresh_bearing_hits"],
            receipt["fresh_anchor_hits"]) == (64, 0, 0)
    assert (receipt["steered_known_scalar_queries"], receipt["steered_raw_hits"],
            receipt["steered_unknown_row_rank"]) == (64, 64, 0)
    lower = int(receipt["ideal_uniform_accepted_hit_bound_denominator"])
    assert receipt["ideal_uniform_accepted_hit_bound_numerator"] == 143
    assert receipt["hit_bound_lt_two_pow_minus_182"] is True
    assert 143 * (1 << 182) < lower
    assert "***" not in gp_log
    match = re.fullmatch(
        r"SPLIT_SUPPORT_GP_PASS triples=(\d+) bearing=(\d+) anchored=(\d+) "
        r"overlap=(\d+) repeated=(\d+) bearing_cells=(\[.*?\]) "
        r"anchored_cells=(\[.*?\]) fresh_bearing=(\d+) fresh_anchor=(\d+)",
        gp_log.strip())
    assert match, gp_log
    assert list(map(int, match.group(1, 2, 3, 4, 5, 8, 9))) == [1540, 571, 969, 0, 0, 0, 0]
    order = ["++", "+-", "-+", "--"]
    for name, text in [("bearing", match.group(6)), ("anchored", match.group(7))]:
        expected = [receipt["character_cells"][name][label] for label in order]
        assert json.loads(text) == expected
    print(json.dumps({"status": "SPLIT_SUPPORT_INDEPENDENT_PASS",
                      "selected_classes": 20, "bearing_targets": 571,
                      "accepted_bearing_targets": 143,
                      "fresh_bearing_hits": 0}, sort_keys=True))


if __name__ == "__main__":
    main()
