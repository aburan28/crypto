"""Replay an N53 scalar certificate in Sage from canonical public inputs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main(candidate_path: Path, workload_path: Path, report_path: Path,
         output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    candidate_raw = candidate_path.read_bytes()
    candidate = json.loads(candidate_raw)
    candidate_id = candidate.pop("candidate_id")
    assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
    workload_raw = workload_path.read_bytes()
    workload = json.loads(workload_raw)
    workload_id = workload["workload_id"]
    record = workload["record"]
    assert workload["identity_sha256"] == sha(canonical(record))
    assert workload_id == workload["identity_sha256"][:12]
    assert candidate["curve"]["curve_id"] == record["curve_id"]
    assert candidate["field"]["degree_n"] == 53
    assert candidate["field"]["characteristic_p"] == 2
    assert candidate["field"]["representation"] == "polynomial_basis"
    report_packed = report_path.read_bytes()
    report_raw = gzip.decompress(report_packed) if report_path.suffix == ".gz" else report_packed
    report = json.loads(report_raw)
    assert report["target"] == record["target_point"]
    assert report["group_verified"] is True
    assert report["recovered_scalar"] is not None

    poly_ring = PolynomialRing(GF(2), "t")
    t = poly_ring.gen()
    modulus = sum((t ** int(power) for power in
                   candidate["field"]["irreducible_polynomial_powers"]),
                  poly_ring.zero())
    field = GF(2**53, "z", modulus=modulus)
    z = field.gen()

    def element(encoded: int):
        assert 0 <= encoded < 2**53
        result = field.zero()
        power = field.one()
        while encoded:
            if encoded & 1:
                result += power
            power *= z
            encoded >>= 1
        return result

    coefficients = candidate["curve"]["coefficients"]
    curve = EllipticCurve(field, [element(coefficients[key]) for key in
                                  ("a1", "a2", "a3", "a4", "a6")])
    generator = curve(*(element(word) for word in candidate["curve"]["generator_G"]))
    target = curve(*(element(word) for word in record["target_point"]))
    subgroup_order = ZZ(candidate["curve"]["subgroup_order_r"])
    recovered = ZZ(report["recovered_scalar"])
    generator_order_ok = subgroup_order * generator == curve(0)
    replay_ok = recovered * generator == target
    assert generator_order_ok and replay_ok
    result = {
        "kind": "s3_n53_independent_sage_scalar_replay_v1",
        "candidate_id": candidate_id,
        "workload_id": workload_id,
        "candidate_manifest_sha256": sha(candidate_raw),
        "workload_manifest_sha256": sha(workload_raw),
        "report_raw_sha256": sha(report_raw),
        "report_gzip_sha256": sha(report_packed) if report_path.suffix == ".gz" else None,
        "subgroup_order_r": int(subgroup_order),
        "recovered_scalar": int(recovered),
        "generator_order_verified": generator_order_ok,
        "scalar_times_generator_equals_public_target": replay_ok,
        "independent_implementation": "Sage finite field and elliptic curve arithmetic",
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--workload", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.candidate, args.workload, args.report, args.output)
