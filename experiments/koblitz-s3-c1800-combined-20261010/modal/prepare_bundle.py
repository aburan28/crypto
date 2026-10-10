"""Package one frozen public N53 point and exact C1200/C1800 expectations."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STUDY = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/CONFIRM_FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def make_bundle() -> dict:
    study = json.loads(STUDY.read_bytes())
    assert study["kind"] == "s3_n53_c1800_combined_confirmation_freeze_v1"
    target = next(row for row in study["targets"] if row["label"] == "T205")
    point = ROOT / target["point_file_relative"]
    assert sha(point.read_bytes()) == target["point_sha256"]
    assert json.loads(point.read_bytes()) == target["public_point"]
    arms = {}
    for local_arm, remote_arm in (("baseline", "baseline"), ("combined", "expanded")):
        frozen = study["arms"][local_arm]
        arms[remote_arm] = {key: frozen[key] for key in (
            "candidate_id", "columns", "prepare_rank", "factor_base_points",
            "factor_base_digest", "root_table_slots", "source_sha256",
            "target_search_schedule", "target_search_prefix_columns")}
    return {
        "kind": "s3_modal_c1800_public_t205_bundle_v1",
        "study_freeze_sha256": sha(STUDY.read_bytes()),
        "relation_seed": study["relation_seed"],
        "relation_workers": study["relation_workers"],
        "target_workers": study["target_workers"],
        "run_order": [[53, "baseline"], [53, "expanded"],
                      [53, "expanded"], [53, "baseline"]],
        "curves": [{
            "n": 53,
            "curve_id": "EC1N53Ckb1h738a572b91f2",
            "label": "T205",
            "public_point": target["public_point"],
            "point_sha256": target["point_sha256"],
            "workload_id": target["workload_id"],
            "arms": arms,
        }],
    }


if __name__ == "__main__":
    output = HERE / "public_bundle.json"
    if output.exists():
        raise FileExistsError(output)
    raw = canonical(make_bundle()) + b"\n"
    assert b"fixture_scalar" not in raw
    output.write_bytes(raw)
    print(json.dumps({"path": str(output), "bytes": len(raw), "sha256": sha(raw)}))
