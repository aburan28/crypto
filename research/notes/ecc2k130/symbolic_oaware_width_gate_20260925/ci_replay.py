#!/usr/bin/env python3
"""Hash-only preflight, then archive-structure replay if a receipt exists."""

import hashlib
import json
from pathlib import Path

from audit import HERE, canonical, preflight


def main() -> None:
    frozen, inp = preflight()
    receipt_path = HERE / "evidence" / "receipt.json"
    if not receipt_path.exists():
        print("FROZEN_HASHES_PASS; OUTCOME_HELD")
        return
    raw = receipt_path.read_bytes()
    assert len(raw) <= inp["caps"]["receipt_bytes"]
    r = json.loads(raw)
    assert raw == canonical(r)
    assert r["schema"] == "symbolic-oaware-width-gate-receipt-v1"
    assert r["frozen_sha256"] == hashlib.sha256((HERE / "FROZEN.json").read_bytes()).hexdigest()
    assert r["source_commit"] == frozen["parent_commit"]
    assert r["toy_mu4_map"]["status"] == "PASS"
    assert r["toy_mu4_map"]["source_affine_points"] == r["toy_mu4_map"]["mu4_affine_chart_points"]
    assert r["toy_mu4_map"]["projective_points"] == r["toy_mu4_map"]["source_affine_points"] + 1
    assert r["toy_mu4_map"]["controls"] == {name: "PASS" for name in inp["toy_controls"]}
    assert r["width"]["status"] == "NOT_ADMITTED"
    assert r["width"]["missing_required_artifacts"] == inp["required_for_n131_admission"]
    assert [(a["factor_bits"], a["affine_chain_floor_bits"]) for a in r["width"]["arms"]] == [(131, 1048), (131, 1179), (130, 1178)]
    assert r["cost"]["wall_seconds"] <= inp["caps"]["wall_seconds"]
    assert r["cost"]["peak_rss_bytes"] <= inp["caps"]["rss_bytes"]
    print("ARCHIVE_STRUCTURE_PASS", hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    main()
