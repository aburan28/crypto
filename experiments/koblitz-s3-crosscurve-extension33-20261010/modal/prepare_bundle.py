"""Package four frozen public points for a source-minimal Modal pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
STUDY = ROOT / "experiments/koblitz-s3-crosscurve-extension33-20261010/FREEZE.json"
INDICES = (25, 28)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def block_order(n: int, index: int) -> tuple[str, ...]:
    return (("reference", "prepared", "prepared", "reference")
            if (index + (n == 53)) % 2 else
            ("prepared", "reference", "reference", "prepared"))


def make_bundle() -> dict:
    parent_raw, study_raw = PARENT.read_bytes(), STUDY.read_bytes()
    parent, study = json.loads(parent_raw), json.loads(study_raw)
    assert parent["kind"] == "s3_rank243_matched_curve_fresh_panel_v1"
    assert study["kind"] == "s3_crosscurve_extension33_freeze_v1"
    assert study["old_freeze_sha256"] == sha(parent_raw)
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    cases, order = [], []
    for index in INDICES:
        label = f"T{index:03}"
        for n in ((41, 53) if index % 2 else (53, 41)):
            target = targets[n, label]
            point_file = PARENT.parent / target["point_file_relative"]
            point_raw = point_file.read_bytes()
            assert sha(point_raw) == target["point_sha256"]
            assert json.loads(point_raw) == target["public_point"]
            curve = parent["curves"][str(n)]
            case = {"n": n, "label": label, "curve_id": curve["curve_id"],
                    "factor_base_points": curve["factor_base_points"],
                    "public_point": target["public_point"],
                    "point_sha256": target["point_sha256"],
                    "workload_id": target["workload_id"],
                    "candidate_ids": {"reference": curve["reference_candidate_id"],
                                      "prepared": curve["rank243_candidate_id"]}}
            cases.append(case)
            # Local panels already use R1/R2 for these candidate/workload IDs.
            repeats = {"reference": 2, "prepared": 2}
            for arm in block_order(n, index):
                repeats[arm] += 1
                order.append({"n": n, "label": label, "arm": arm,
                              "repetition": repeats[arm],
                              "candidate_id": case["candidate_ids"][arm],
                              "workload_id": case["workload_id"],
                              "run_id": f"{case['candidate_ids'][arm]}W{case['workload_id']}R{repeats[arm]}"})
    assert len(cases) == 4 and len(order) == 16
    return {"kind": "s3_crosscurve_modal_public_pilot_bundle_v1",
            "parent_freeze_sha256": sha(parent_raw),
            "extension_freeze_sha256": sha(study_raw),
            "source_commit": study["source_commit"],
            "relation_seed": parent["relation_seed"], "workers": parent["workers"],
            "prepare_rank": 243, "orbit_columns": 244,
            "cases": cases, "run_order": order}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    raw = canonical(make_bundle()) + b"\n"
    assert b"fixture_scalar" not in raw and b"scalar_derivation" not in raw
    args.output.write_bytes(raw)
    print(json.dumps({"bytes": len(raw), "sha256": sha(raw),
                      "cases": 4, "runs": 16}, sort_keys=True))
