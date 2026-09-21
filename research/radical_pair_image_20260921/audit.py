#!/usr/bin/env python3
"""Verify frozen radical pair-image evidence and reporting boundaries."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    args = parser.parse_args()
    run = Path(args.run).resolve()
    root = Path(__file__).resolve().parent
    manifest = json.loads((run / "manifest.json").read_text())
    processes = [json.loads(line) for line in (run / "processes.jsonl").read_text().splitlines()]
    certificates = json.loads((run / "certificates.json").read_text())
    cells = json.loads((run / "cells.json").read_text())
    summary = json.loads((run / "summary.json").read_text())

    source_checks = {}
    for name, expected in manifest["source_hashes"].items():
        actual = sha256_file(root / name)
        source_checks[name] = {"expected": expected, "actual": actual, "match": actual == expected}
    assert all(item["match"] for item in source_checks.values())

    ordinals = [row["ordinal"] for row in processes]
    assert ordinals == list(range(len(processes)))
    assert len(processes) == manifest["processes"]
    assert all(row["status"] == "COMPLETED" and row["returncode"] == 0
               for row in processes)
    assert all(row["max_f4_degree"] is not None and row["rounds"]
               for row in processes)
    for row in processes:
        assert sha256_file(run / row["input"]) == row["input_sha256"]
        assert sha256_file(run / row["stdout"]) == row["stdout_sha256"]
        assert sha256_file(run / row["stderr"]) == row["stderr_sha256"]
        assert row["certificate_key"] in certificates

    assert all(cert["presentation_equivalence"] and cert["nondegenerate_group_equivalence"]
               for cert in certificates.values())
    for key, cert in certificates.items():
        if key.startswith("binary-"):
            assert cert["unordered_pair_enumeration"]
            assert cert["quadratic_root_recovery"]
            assert cert["polynomial_evaluation"]
            assert cert["root_recovery"]
    assert min(cert["full_group_coverage"] for cert in certificates.values()) >= 0.8
    assert summary["all_processes_completed"] and summary["all_certificates_pass"]
    assert summary["primary_metrics"]["total_common_operations"] is None
    assert summary["primary_metrics"]["S"] is None
    assert summary["primary_metrics"]["rho_ratio"] is None
    assert summary["primary_metrics"]["generic_floor_ratio"] is None

    algebraic_boundaries = []
    for cell in cells:
        if cell["label"].startswith("prime"):
            b = cell["factor_base_size"]
            preprocessing = cell["radical_preprocessing"]
            assert preprocessing["basis_size"] == b + 1
            assert preprocessing["basis_degree"] == b
            algebraic_boundaries.append({
                "cell": cell["label"],
                "b": b,
                "ordered_membership_regularity": 2 * b - 1,
                "radical_membership_regularity": b,
                "direct_over_floor": (2 * b - 1) / b,
                "radical_over_floor": 1.0,
            })

    final = {
        "status": "pass",
        "processes": len(processes),
        "cells": len(cells),
        "certificates": len(certificates),
        "source_checks": source_checks,
        "algebraic_boundaries": algebraic_boundaries,
        "claim_boundary": "solver-stage engineering diagnostic only",
        "end_to_end_metrics_null": True,
    }
    (run / "final-audit.json").write_text(json.dumps(final, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "processes": len(processes),
                      "certificates": len(certificates)}, indent=2))


if __name__ == "__main__":
    main()
