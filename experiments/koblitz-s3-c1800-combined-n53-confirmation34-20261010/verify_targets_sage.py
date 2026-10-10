"""Independently audit the frozen N53 one-target fixtures with checked Sage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GENERATED = HERE / "GENERATED_TARGETS.json"
RUNTIME = HERE / "sage-runtime-info.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def verify(output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    generated_raw = GENERATED.read_bytes()
    generated = json.loads(generated_raw)
    parent_raw = PARENT.read_bytes()
    parent = json.loads(parent_raw)
    runtime_raw = RUNTIME.read_bytes()
    assert json.loads(runtime_raw)["status"] == "verified"
    assert generated["kind"] == "s3_n53_c1800_combined_confirmation34_public_targets_v1"
    assert generated["parent_freeze_sha256"] == sha(parent_raw)
    assert generated["fixture_generator_sha256"] == sha((HERE / "make_targets_sage.py").read_bytes())
    assert generated["sage_runtime_info_sha256"] == sha(runtime_raw)
    targets = generated["targets"]
    assert [item["label"] for item in targets] == [f"T{i:03}" for i in range(247, 281)]
    assert len({item["workload_id"] for item in targets}) == len(targets) == 34

    declared = parent["curves"]["53"]
    ring = PolynomialRing(GF(2), "t")
    t = ring.gen()
    field = GF(2**53, "z", modulus=sum(
        (t ** int(power) for power in declared["modulus_powers"]), ring.zero()))
    z = field.gen()

    def decode(word: int):
        assert 0 <= word < 2**53
        return sum((z**i for i in range(53) if (word >> i) & 1), field.zero())

    curve = EllipticCurve(field, [decode(word) for word in (1, 0, 0, 0, 1)])
    generator = curve(*(decode(word) for word in declared["generator"]))
    r = ZZ(declared["r"])
    assert r * generator == curve(0)
    seen_points = set()
    rows = []
    for item in targets:
        index = int(item["label"][1:])
        scalar_input = (b"s3-pair-target-v1:" + parent["seed_label"].encode() + b":" +
                        struct.pack(">II", 53, index))
        derived = sha(scalar_input)
        scalar = 1 + int(derived, 16) % (int(r) - 1)
        assert item["scalar_derivation_sha256"] == derived
        assert item["fixture_scalar"] == scalar
        point = item["public_point"]
        assert isinstance(point, list) and len(point) == 2
        assert tuple(point) not in seen_points
        seen_points.add(tuple(point))
        point_file = ROOT / item["point_file_relative"]
        point_raw = point_file.read_bytes()
        assert sha(point_raw) == item["point_sha256"]
        assert json.loads(point_raw) == point
        q = curve(decode(point[0]), decode(point[1]))
        assert scalar * generator == q
        assert r * q == curve(0)

        workload_file = ROOT / item["workload_manifest_relative"]
        workload_raw = workload_file.read_bytes()
        assert sha(workload_raw) == item["workload_manifest_sha256"]
        workload = json.loads(workload_raw)
        assert workload["record"]["curve_id"] == declared["curve_id"]
        assert workload["record"]["target_point"] == point
        assert workload["record"]["target_count"] == 1
        assert workload["record"]["target_index"] == index
        assert workload["identity_sha256"] == sha(canonical(workload["record"]))
        assert workload["workload_id"] == workload["identity_sha256"][:12]
        assert workload["workload_id"] == item["workload_id"]
        assert workload["identity_sha256"] == item["workload_identity_sha256"]
        rows.append({"label": item["label"], "workload_id": item["workload_id"],
                     "point_sha256": item["point_sha256"], "scalar_replay": True})

    receipt = {
        "kind": "s3_n53_c1800_combined_confirmation34_target_verification_v1",
        "generated_targets_sha256": sha(generated_raw),
        "parent_freeze_sha256": sha(parent_raw),
        "runtime_info_sha256": sha(runtime_raw),
        "verifier_sha256": sha(Path(__file__).read_bytes()),
        "target_count": len(rows),
        "all_scalar_replays_verified": True,
        "rows": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"target_count": len(rows), "all_scalar_replays_verified": True,
                      "generated_targets_sha256": sha(generated_raw)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    verify(parser.parse_args().output.resolve())
