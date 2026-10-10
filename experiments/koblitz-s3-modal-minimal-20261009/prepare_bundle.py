"""Freeze the public N41/N53 T181 points and exact executable expectations."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
STUDIES = {
    41: ROOT / "experiments/koblitz-s3-n41-c800-fullrank-20261009/FREEZE.json",
    53: ROOT / "experiments/koblitz-s3-target-conjugate-20261009/C1200_FULLRANK_FREEZE.json",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def make_bundle() -> dict:
    parent = json.loads(PARENT.read_bytes())
    curves = []
    for n in (41, 53):
        study = json.loads(STUDIES[n].read_bytes())
        target = next(item for item in parent["targets"]
                      if item["n"] == n and item["label"] == "T181")
        point_path = PARENT.parent / target["point_file_relative"]
        point_raw = point_path.read_bytes()
        assert sha(point_raw) == target["point_sha256"]
        assert json.loads(point_raw) == target["public_point"]
        assert study["curve_id"] == parent["curves"][str(n)]["curve_id"]
        arms = {}
        for arm in ("baseline", "expanded"):
            frozen = study["arms"][arm]
            assert frozen["source_sha256"] in (
                "71f527cf7d4f401ad252a6f5ad3fa504c3f8becac9567e254e10d0b48cd75cfa",
                "0ce2fa5f52c15776fa9cd80c173e14c3c405815b6fb2628d6e75d304b7af2559",
            )
            arms[arm] = {
                "candidate_id": frozen["candidate_id"],
                "columns": frozen["columns"],
                "prepare_rank": frozen["requested_pre_target_rank"],
                "factor_base_points": frozen["factor_base_points"],
                "factor_base_digest": frozen["factor_base_digest"],
                "root_table_slots": frozen["root_table_slots"],
                "source_sha256": frozen["source_sha256"],
            }
        curves.append({
            "n": n, "curve_id": study["curve_id"], "label": "T181",
            "public_point": target["public_point"],
            "point_sha256": target["point_sha256"],
            "workload_id": target["workload_id"],
            "arms": arms,
        })
    return {
        "kind": "s3_modal_minimal_public_t181_bundle_v1",
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "study_freeze_sha256": {str(n): sha(STUDIES[n].read_bytes()) for n in (41, 53)},
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "run_order": [[41, "baseline"], [41, "expanded"],
                      [41, "expanded"], [41, "baseline"],
                      [53, "expanded"], [53, "baseline"],
                      [53, "baseline"], [53, "expanded"]],
        "curves": curves,
    }


if __name__ == "__main__":
    output = HERE / "public_bundle.json"
    if output.exists():
        raise FileExistsError(output)
    bundle = make_bundle()
    raw = canonical(bundle) + b"\n"
    assert b"fixture_scalar" not in raw
    output.write_bytes(raw)
    print(json.dumps({"path": str(output), "bytes": len(raw), "sha256": sha(raw)}))
