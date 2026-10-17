#!/usr/bin/env python3
"""Check input custody and exact agreement with the PARI support replay."""

import collections
import hashlib
import json
import re
import sys
from pathlib import Path


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def coeff(values):
    return [str(value[0] if isinstance(value, list) else value) for value in values]


def main():
    if len(sys.argv) != 8:
        raise SystemExit("usage: verify_divisor_support_receipt.py POS NEG GUIDED_NEG RECEIPT GP_INPUT GP_LOG PRODUCER")
    raw = [Path(name).read_bytes() for name in sys.argv[1:4]]
    docs = list(map(json.loads, raw))
    receipt = json.loads(Path(sys.argv[4]).read_bytes())
    gp_input = Path(sys.argv[5]).read_bytes()
    gp_log = Path(sys.argv[6]).read_text()
    assert receipt["input_sha256"] == list(map(sha, raw))
    assert receipt["gp_input_sha256"] == sha(gp_input)
    assert receipt["source_sha256"] == sha(Path(sys.argv[7]).read_bytes())
    assert receipt["status"] == "exact_fixed_origin_divisor_support_audit"
    assert (receipt["class_count"], receipt["squarefree_cubics"]) == (160, 160)
    assert (receipt["distinct_numerator_cubic_sections"],
            receipt["distinct_denominator_cubic_sections"]) == (160, 160)
    assert receipt["pairwise_x_gcd_degree_counts"] == {"0": 12720}
    assert receipt["first_shared_x_witnesses"] == []
    positive, negative, guided_negative = [doc["representatives"] for doc in docs]
    groups = [
        ("base_positive", positive[:16]),
        ("guided_positive", positive[16:]),
        ("base_negative", negative),
        ("guided_negative", guided_negative),
    ]
    assert [len(rows) for _, rows in groups] == [16, 64, 16, 64]
    records = receipt["records"]
    assert len(records) == 160
    offsets = 0
    expected_counts = []
    for group, rows in groups:
        counts = collections.Counter()
        for index, source in enumerate(rows):
            record = records[offsets + index]
            assert record["group"] == group and record["index"] == index
            assert record["x_coefficients"] == coeff(source["representative_x_polynomial"])
            assert record["y_coefficients"] == coeff(source["representative_y_polynomial"])
            counts[record["linear_roots"]] += 1
        assert receipt["linear_root_count_by_group"][group] == {
            str(key): value for key, value in sorted(counts.items())}
        expected_counts.append([counts.get(i, 0) for i in range(4)])
        offsets += len(rows)
    assert "***" not in gp_log
    line = gp_log.strip()
    match = re.fullmatch(
        r"DIVISOR_SUPPORT_GP_PASS classes=(\d+) squarefree=(\d+) "
        r"pairs=(\d+) shared_x_pairs=(\d+) root_counts=(\[.*\])", line)
    assert match, line
    assert list(map(int, match.groups()[:4])) == [160, 160, 12720, 0]
    assert json.loads(match.group(5)) == expected_counts
    print(json.dumps({"status": "DIVISOR_SUPPORT_INDEPENDENT_PASS",
                      "classes": 160, "disjoint_pairs": 12720,
                      "split_cubics": sum(group[3] for group in expected_counts),
                      "root_counts": expected_counts}, sort_keys=True))


if __name__ == "__main__":
    main()
