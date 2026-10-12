#!/usr/bin/env python3
"""Freeze the preregistered scalar stream and #737 point-only holdouts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PR737 = HERE.parent / "paired_fullrank_clean_evidence_20260925/clean_summary.json"
DOMAIN = b"ECC2K-COMPACT-BASE-SWEEP-20260925-v1/"
ORDERS = {37: 230603167, 41: 549756390943}
ARMS = {37: [(1, 146), (2, 1317), (3, 5125)],
        41: [(4, 7), (8, 71), (12, 255)]}
COUNT = 4096


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_once(path: Path, data: bytes) -> None:
    if path.exists():
        assert path.read_bytes() == data, f"frozen input changed: {path}"
    else:
        path.write_bytes(data)


def scalars(n: int, order: int) -> list[int]:
    values, seen = [], set()
    counter = 0
    while len(values) < COUNT:
        digest = hashlib.sha256(DOMAIN + str(n).encode() + b"/"
                                + str(counter).encode()).digest()
        value = 1 + int.from_bytes(digest, "big") % (order - 1)
        counter += 1
        if value not in seen:
            values.append(value)
            seen.add(value)
    return values


def main() -> None:
    summary = json.loads(PR737.read_bytes())
    assert summary["classification"] == "clean main-derived source; paired full-rank control"
    spec = {"schema_version": "1.0", "domain": DOMAIN.decode(),
            "scalar_count": COUNT, "arms": {}, "holdouts": {},
            "pr737_clean_summary_sha256": sha(PR737.read_bytes())}
    for n, order in ORDERS.items():
        stream = scalars(n, order)
        scalar_data = ("".join(f"{value}\n" for value in stream)).encode()
        scalar_path = HERE / f"training_scalars_n{n}.txt"
        write_once(scalar_path, scalar_data)
        rows = [run for run in summary["runs"] if run["n"] == n]
        assert len(rows) == 3
        assert [run["seed"] for run in rows] == [202609250000+n, 202609250100+n, 202609250200+n]
        assert all(run["recovered_scalar"] not in stream for run in rows)
        holdouts = [{"seed": run["seed"], "q": run["q"],
                     "scalar_validator_only": run["recovered_scalar"]}
                    for run in rows]
        point_data = ("".join(json.dumps(run["q"], separators=(",", ":")) + "\n"
                              for run in holdouts)).encode()
        point_path = HERE / f"holdout_points_n{n}.jsonl"
        write_once(point_path, point_data)
        spec["arms"][str(n)] = [{"R": r, "eta_numerator": eta,
                                  "eta_denominator": 1000000,
                                  "F": 2*n*r, "subgroup_order": order}
                                 for r, eta in ARMS[n]]
        spec["holdouts"][str(n)] = holdouts
        spec[f"training_scalars_n{n}_sha256"] = sha(scalar_data)
        spec[f"holdout_points_n{n}_sha256"] = sha(point_data)
    spec_data = (json.dumps(spec, indent=2, sort_keys=True) + "\n").encode()
    write_once(HERE / "input_spec.json", spec_data)
    print(json.dumps({"input_spec_sha256": sha(spec_data),
                      "training_sha256": {n: spec[f"training_scalars_n{n}_sha256"] for n in ORDERS},
                      "holdout_sha256": {n: spec[f"holdout_points_n{n}_sha256"] for n in ORDERS}},
                     sort_keys=True))


if __name__ == "__main__":
    main()
