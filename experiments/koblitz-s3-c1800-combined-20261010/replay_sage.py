"""Independently replay every N53 C1800 combined-kernel target scalar through checked Sage."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLAN = HERE / "CONFIRM_FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def replay(run_dir: Path, runtime_info: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    runtime_raw = runtime_info.read_bytes()
    assert runtime_raw
    plan_raw = PLAN.read_bytes()
    plan = json.loads(plan_raw)
    manifests = {}
    for arm in ("baseline", "combined"):
        frozen = plan["arms"][arm]
        raw = (ROOT / frozen["candidate_manifest_relative"]).read_bytes()
        assert sha(raw) == frozen["candidate_manifest_sha256"]
        manifest = json.loads(raw)
        candidate_id = manifest.pop("candidate_id")
        assert candidate_id == frozen["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(manifest))[:12])
        manifests[arm] = manifest
    assert manifests["baseline"]["field"] == manifests["combined"]["field"]
    assert manifests["baseline"]["curve"] == manifests["combined"]["curve"]
    field_record = manifests["baseline"]["field"]
    curve_record = manifests["baseline"]["curve"]
    poly_ring = PolynomialRing(GF(2), "t")
    t = poly_ring.gen()
    modulus = sum((t ** int(power) for power in
                   field_record["irreducible_polynomial_powers"]), poly_ring.zero())
    field = GF(2 ** field_record["degree_n"], "z", modulus=modulus)
    z = field.gen()

    def element(encoded: int):
        assert 0 <= encoded < 2 ** field_record["degree_n"]
        result, power = field.zero(), field.one()
        while encoded:
            if encoded & 1:
                result += power
            power *= z
            encoded >>= 1
        return result

    coefficients = curve_record["coefficients"]
    curve = EllipticCurve(field, [element(coefficients[key]) for key in
                                  ("a1", "a2", "a3", "a4", "a6")])
    generator = curve(*(element(word) for word in curve_record["generator_G"]))
    subgroup_order = ZZ(curve_record["subgroup_order_r"])
    assert subgroup_order * generator == curve(0)

    targets = {target["label"]: target for target in plan["targets"]}
    rows = []
    for line in (run_dir / "ledger.jsonl").read_bytes().splitlines():
        entry = json.loads(line)
        target = targets[entry["label"]]
        receipt_raw = (run_dir / entry["relative_receipt"]).read_bytes()
        assert sha(receipt_raw) == entry["receipt_sha256"]
        receipt = json.loads(receipt_raw)
        assert receipt["verified_native_and_fixture"]
        assert receipt["candidate_id"] == plan["arms"][entry["arm"]]["candidate_id"]
        assert receipt["workload_id"] == target["workload_id"]
        report_packed = (run_dir / entry["relative_receipt"]).parent / "report.json.gz"
        packed = report_packed.read_bytes()
        report_raw = gzip.decompress(packed)
        assert sha(packed) == receipt["archive"]["report.json"]["gzip_sha256"]
        report = json.loads(report_raw)
        assert report["target"] == target["public_point"]
        scalar = ZZ(report["recovered_scalar"])
        assert scalar == target["fixture_scalar"]
        public_point = curve(*(element(word) for word in target["public_point"]))
        assert scalar * generator == public_point
        rows.append({"label": entry["label"], "arm": entry["arm"],
                     "candidate_id": receipt["candidate_id"],
                     "workload_id": receipt["workload_id"],
                     "report_raw_sha256": sha(report_raw),
                     "recovered_scalar": int(scalar),
                     "scalar_times_generator_equals_public_target": True})
    assert len(rows) == len(plan["target_labels"]) * 4
    output.write_text(json.dumps({
        "kind": "s3_n53_c1800_combined_confirmation_sage_replay_v1",
        "freeze_sha256": sha(plan_raw),
        "runtime_info_sha256": sha(runtime_raw),
        "generator_order_verified": True,
        "verified_run_count": len(rows),
        "runs": rows,
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verified_run_count": len(rows),
                      "generator_order_verified": True}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--runtime-info", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    replay(args.run_dir, args.runtime_info, args.output)
