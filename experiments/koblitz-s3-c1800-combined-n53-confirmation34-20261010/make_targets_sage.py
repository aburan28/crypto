"""Freeze 34 fresh N53 public points T247-T280 with checked Sage."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import struct

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PREVIOUS = ROOT / "experiments/koblitz-s3-c1800-confirmation34-20261010/GENERATED_TARGETS.json"
COMBINED = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/GENERATED_TARGETS.json"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    parent_raw = PARENT.read_bytes()
    parent = json.loads(parent_raw)
    assert parent["input_law"] == (
        "1 + int_be(SHA256('s3-pair-target-v1:' || seed_label || ':' || "
        "u32be(n) || u32be(index))) mod (r-1)")
    declared = parent["curves"]["53"]
    r = declared["r"]
    ring = PolynomialRing(GF(2), "t")
    t = ring.gen()
    modulus = sum((t ** int(k) for k in declared["modulus_powers"]), ring.zero())
    field = GF(2**53, "z", modulus=modulus)
    z = field.gen()

    def decode(word: int):
        value, power = field.zero(), field.one()
        while word:
            if word & 1:
                value += power
            power *= z
            word >>= 1
        return value

    def encode(value) -> int:
        return sum(int(bit) << index for index, bit in
                   enumerate(value.polynomial().list()))

    curve = EllipticCurve(field, [decode(word) for word in (1, 0, 0, 0, 1)])
    generator = curve(*(decode(word) for word in declared["generator"]))
    assert ZZ(r) * generator == curve(0)

    def fixture(index: int) -> tuple[int, str, list[int]]:
        raw = (b"s3-pair-target-v1:" + parent["seed_label"].encode() + b":" +
               struct.pack(">II", 53, index))
        derived = sha(raw)
        scalar = 1 + int(derived, 16) % (r - 1)
        point = ZZ(scalar) * generator
        return scalar, derived, [encode(point[0]), encode(point[1])]

    existing = next(row for row in parent["targets"]
                    if row["n"] == 53 and row["label"] == "T200")
    assert fixture(200) == (existing["fixture_scalar"],
                            existing["scalar_derivation_sha256"],
                            existing["public_point"])
    previous_raw = PREVIOUS.read_bytes()
    previous = json.loads(previous_raw)
    combined_raw = COMBINED.read_bytes()
    combined = json.loads(combined_raw)
    assert [row["label"] for row in previous["targets"]] == [
        f"T{index:03}" for index in range(213, 247)]
    assert [row["label"] for row in combined["targets"]] == [
        f"T{index:03}" for index in range(205, 213)]
    prior = [row for row in parent["targets"] if row["n"] == 53]
    prior.extend(combined["targets"])
    prior.extend(previous["targets"])
    assert fixture(246) == (previous["targets"][-1]["fixture_scalar"],
                            previous["targets"][-1]["scalar_derivation_sha256"],
                            previous["targets"][-1]["public_point"])
    seen_points = {tuple(row["public_point"]) for row in prior}
    seen_workloads = {row["workload_id"] for row in prior
                      if "workload_id" in row}
    assert len(seen_points) == len(prior)

    records = []
    for index in range(247, 281):
        label = f"T{index:03}"
        scalar, derived, public_point = fixture(index)
        assert tuple(public_point) not in seen_points
        seen_points.add(tuple(public_point))
        point_rel = f"fixtures/{label}/public_target.json"
        point_path = HERE / point_rel
        point_path.parent.mkdir(parents=True, exist_ok=True)
        point_raw = json.dumps(public_point, separators=(",", ":")).encode() + b"\n"
        if point_path.exists():
            assert point_path.read_bytes() == point_raw
        else:
            point_path.write_bytes(point_raw)
        workload_record = {
            "curve_id": declared["curve_id"],
            "preparation_state_at_online_start": "candidate_policy_reusable_setup_ready",
            "relation_seed": parent["relation_seed"],
            "scalar_derivation_sha256": derived,
            "subgroup_order_r": r,
            "target_count": 1,
            "target_index": index,
            "target_input_law": "known_scalar_sha256_public_point_fixture_v1",
            "target_label": label,
            "target_point": public_point,
            "target_seed_label": parent["seed_label"],
        }
        identity = sha(canonical(workload_record))
        workload_id = identity[:12]
        assert workload_id not in seen_workloads
        seen_workloads.add(workload_id)
        manifest = {"identity_sha256": identity, "record": workload_record,
                    "workload_id": workload_id}
        workload_path = HERE / "workloads" / f"W{workload_id}.json"
        workload_path.parent.mkdir(parents=True, exist_ok=True)
        workload_raw = json.dumps(manifest, sort_keys=True, indent=2).encode() + b"\n"
        if workload_path.exists():
            assert workload_path.read_bytes() == workload_raw
        else:
            workload_path.write_bytes(workload_raw)
        records.append({"label": label, "n": 53, "fixture_scalar": scalar,
                        "scalar_derivation_sha256": derived,
                        "public_point": public_point,
                        "point_file_relative": str(point_path.relative_to(ROOT)),
                        "point_sha256": sha(point_raw),
                        "workload_manifest_relative": str(workload_path.relative_to(ROOT)),
                        "workload_manifest_sha256": sha(workload_raw),
                        "workload_id": workload_id,
                        "workload_identity_sha256": identity})
    runtime_info = HERE / "sage-runtime-info.json"
    assert runtime_info.exists()
    output = HERE / "GENERATED_TARGETS.json"
    output_raw = (json.dumps({
        "kind": "s3_n53_c1800_combined_confirmation34_public_targets_v1",
        "parent_freeze_sha256": sha(parent_raw),
        "previous_targets_sha256": sha(previous_raw),
        "combined_targets_sha256": sha(combined_raw),
        "fixture_generator_sha256": sha(Path(__file__).read_bytes()),
        "sage_runtime_info_sha256": sha(runtime_info.read_bytes()),
        "generator_order_verified": True,
        "existing_t200_fixture_replayed": True,
        "existing_t246_fixture_replayed": True,
        "targets": records,
    }, sort_keys=True, indent=2) + "\n").encode()
    if output.exists():
        assert output.read_bytes() == output_raw
    else:
        output.write_bytes(output_raw)
    print(json.dumps({"labels": [row["label"] for row in records],
                      "workload_ids": [row["workload_id"] for row in records]}))


if __name__ == "__main__":
    main()
