"""Independently replay the four N53 striped-schedule scalars in Sage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "koblitz-s3-c1800-matchedsource-n53-20261010"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main(stripe_dir: Path, runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_bytes())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    diagnostic_path = HERE / "CALL_DIAGNOSTIC.json"
    diagnostic = json.loads(diagnostic_path.read_bytes())
    parent_raw = (PARENT / "FREEZE.json").read_bytes()
    assert diagnostic["parent_freeze_sha256"] == sha(parent_raw)
    parent = json.loads(parent_raw)
    targets = {target["label"]: target for target in parent["selected_targets"]}
    manifest_path = HERE / diagnostic["candidate_manifest_relative"]
    assert sha(manifest_path.read_bytes()) == diagnostic["candidate_manifest_sha256"]
    candidate = json.loads(manifest_path.read_bytes())
    assert candidate["candidate_id"] == diagnostic["candidate_id"]
    curve_record = candidate["curve"]
    assert curve_record["subgroup_order_r"] == 21044858204113
    assert curve_record["cofactor"] == 428
    ring = PolynomialRing(GF(2), "t")
    t = ring.gen()
    field = GF(2**53, name="z", modulus=t**53 + t**6 + t**2 + t + 1)
    powers = [field.gen()**bit for bit in range(53)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << 53
        value = field(0)
        while word:
            bit = (word & -word).bit_length() - 1
            value += powers[bit]
            word &= word - 1
        return value

    curve = EllipticCurve(field, [1, 0, 0, 0, 1])

    def point(words: list[int]):
        assert len(words) == 2
        return curve(from_word(words[0]), from_word(words[1]))

    generator = point(curve_record["generator_G"])
    order = curve_record["subgroup_order_r"]
    assert generator != curve(0) and order * generator == curve(0)
    checked = []
    for row in diagnostic["rows"]:
        label = row["label"]
        target = targets[label]
        raw = (stripe_dir / f"s3-c1800-stripe-smoke-{label}.json").read_bytes()
        report = json.loads(raw)
        assert sha(raw) == row["report_sha256"]
        assert row["workload_id"] == target["workload_id"]
        assert row["point_sha256"] == target["point_sha256"]
        assert report["target"] == target["public_point"]
        assert report["factor_base_digest"] == candidate["factor_base"]["enumerated_set_digest"]
        assert report["group_verified"] is True
        scalar = int(report["recovered_scalar"])
        assert 0 < scalar < order
        assert scalar == row["recovered_scalar"] == target["fixture_scalar"]
        assert scalar * generator == point(target["public_point"])
        checked.append({"label": label, "run_id": row["run_id"],
                        "report_sha256": sha(raw), "recovered_scalar": scalar,
                        "sage_scalar_replay": True})
    assert len(checked) == 4
    receipt = {
        "kind": "s3_c1800_stripe_call_diagnostic_sage_replay_v1",
        "diagnostic_sha256": sha(diagnostic_path.read_bytes()),
        "replay_script_sha256": sha(Path(__file__).read_bytes()),
        "sage_runtime_info_sha256": sha(runtime_path.read_bytes()),
        "sage_version": SAGE_VERSION,
        "run_count": len(checked), "all_scalars_independently_replayed": True,
        "runs": checked,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"runs": len(checked), "independent_replay": True}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stripe-dir", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.stripe_dir.resolve(), args.runtime_info.resolve(), args.receipt.resolve())
