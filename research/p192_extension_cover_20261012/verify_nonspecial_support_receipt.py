#!/usr/bin/env python3
"""Bind the independent PARI noncollinearity replay to the frozen support."""

import hashlib
import json
import sys
from pathlib import Path


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    if len(sys.argv) != 7:
        raise SystemExit("usage: verify_nonspecial_support_receipt.py SUPPORT RECEIPT GP_INPUT GP_LOG PRODUCER GP_SOURCE")
    support_path, receipt_path, input_path, log_path, producer_path, gp_source_path = map(Path, sys.argv[1:])
    support_raw = support_path.read_bytes()
    support = json.loads(support_raw)
    receipt = json.loads(receipt_path.read_bytes())
    gp_raw = input_path.read_bytes()
    gp_log = log_path.read_text()
    assert support["status"] == "exact_fixed_origin_divisor_support_audit"
    assert (support["class_count"], support["squarefree_cubics"]) == (160, 160)
    assert receipt["status"] == "fixed_origin_nonspecial_coefficient_screen"
    assert receipt["support_sha256"] == sha(support_raw)
    assert receipt["source_sha256"] == sha(producer_path.read_bytes())
    assert receipt["gp_input_sha256"] == sha(gp_raw)
    assert receipt["field_prime"] == support["field_prime"]
    assert receipt["source_label"] == support_path.name
    assert receipt["class_count"] == receipt["squarefree_distinct_abscissa_classes"] == 160
    assert receipt["quadratic_y_coefficient_nonzero"] == 160
    assert receipt["collinearity_witness_indices"] == []
    polys = []
    for record in support["records"]:
        values = record["y_coefficients"]
        assert len(values) == 3 and int(values[2]) != 0
        polys.append(" + ".join(f"{int(c)}*z^{i}" for i, c in enumerate(values)))
    expected_gp = (f"P={support['field_prime']};\nY=[" + ",".join(polys) + "];\n").encode()
    assert gp_raw == expected_gp
    assert "NONSPECIAL_SUPPORT_GP_PASS classes=160 nonzero_quadratic=160" == gp_log.strip()
    assert "***" not in gp_log
    assert b"NONSPECIAL_SUPPORT_GP_PASS" in gp_source_path.read_bytes()
    print("NONSPECIAL_SUPPORT_INDEPENDENT_PASS classes=160")


if __name__ == "__main__":
    main()
