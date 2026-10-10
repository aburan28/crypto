"""Independently verify the C3200 prefetch T025 scalar and subgroup."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FREEZE = HERE / "CONTROL_FREEZE.json"
BASE_MANIFEST = ROOT / "experiments/koblitz-s3-c3200-projective-20261010/identities/candidates/IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0h0e595bd6143a.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def replay(receipt_path: Path, report_path: Path, runtime_path: Path,
           output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    freeze_raw = FREEZE.read_bytes()
    freeze = json.loads(freeze_raw)
    receipt_raw = receipt_path.read_bytes()
    receipt = json.loads(receipt_raw)
    packed = report_path.read_bytes()
    raw = gzip.decompress(packed)
    report = json.loads(raw)
    runtime_raw = runtime_path.read_bytes()
    assert json.loads(runtime_raw)["status"] == "verified"
    assert receipt["freeze_sha256"] == sha(freeze_raw)
    assert receipt["scalar_verified"] and receipt["full_rank_before_target"]
    assert receipt["index_layout_ok"] and receipt["memory_within_envelope"]
    assert receipt["raw_report_sha256"] == sha(raw)
    assert report["target"] == freeze["public_point"]
    assert report["recovered_scalar"] == freeze["fixture_scalar"]
    manifest_raw = BASE_MANIFEST.read_bytes()
    manifest = json.loads(manifest_raw)
    field_record, curve_record = manifest["field"], manifest["curve"]
    polynomial_ring = PolynomialRing(GF(2), "t")
    t = polynomial_ring.gen()
    modulus = sum((t ** int(power) for power in
                   field_record["irreducible_polynomial_powers"]),
                  polynomial_ring.zero())
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
    target = curve(*(element(word) for word in freeze["public_point"]))
    scalar = ZZ(report["recovered_scalar"])
    assert scalar * generator == target
    result = {
        "kind": "s3_n53_c3200_prefetch_control_sage_replay_v1",
        "freeze_sha256": sha(freeze_raw),
        "receipt_sha256": sha(receipt_raw),
        "report_raw_sha256": sha(raw),
        "runtime_info_sha256": sha(runtime_raw),
        "replay_source_sha256": sha(Path(__file__).read_bytes()),
        "generator_order_verified": True,
        "scalar_times_generator_equals_public_target": True,
        "recovered_scalar": int(scalar),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("generator_order_verified",
                       "scalar_times_generator_equals_public_target",
                       "recovered_scalar")}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--runtime-info", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    replay(args.receipt.resolve(), args.report.resolve(),
           args.runtime_info.resolve(), args.output.resolve())
