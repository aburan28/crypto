"""Independently replay all 48 N41 confirmation scalars in checked Sage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main(ledger_path: Path, runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_bytes())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    rows_raw = ledger_path.read_bytes()
    rows = [json.loads(line) for line in rows_raw.splitlines()]
    assert len(rows) == freeze["planned_run_count"] == 48
    assert [(row["label"], row["arm"], row["repetition"]) for row in rows] == [
        (entry["label"], entry["arm"], entry["repetition"])
        for entry in freeze["run_order"]]
    target_map = {target["label"]: target for target in freeze["selected_targets"]}
    curve_record = json.loads((HERE / freeze["arms"]["baseline"]["candidate_manifest_relative"]).read_bytes())["curve"]
    assert curve_record["subgroup_order_r"] == 549756390943
    assert curve_record["cofactor"] == 4
    ring = PolynomialRing(GF(2), "t")
    t = ring.gen()
    field = GF(2**41, name="z", modulus=t**41 + t**3 + 1)
    powers = [field.gen()**bit for bit in range(41)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << 41
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
    subgroup_order = curve_record["subgroup_order_r"]
    assert generator != curve(0) and subgroup_order * generator == curve(0)
    checked = []
    for row in rows:
        target = target_map[row["label"]]
        arm = freeze["arms"][row["arm"]]
        stem = f"{row['label']}_{row['arm']}_R{row['repetition']}"
        report_raw = (ledger_path.parent / f"{stem}.json").read_bytes()
        report = json.loads(report_raw)
        assert row["exit_code"] == 0 and row["verified_native_and_fixture"] is True
        assert row["report_sha256"] == sha(report_raw)
        assert row["candidate_id"] == arm["candidate_id"]
        assert row["workload_id"] == target["workload_id"]
        assert report["group_verified"] is True
        assert report["target"] == target["public_point"]
        assert report["factor_base_digest"] == arm["factor_base_digest"]
        assert len(report["target_relation_indices"]) == 4
        scalar = int(report["recovered_scalar"])
        assert 0 < scalar < subgroup_order
        assert scalar == target["fixture_scalar"]
        assert scalar * generator == point(target["public_point"])
        checked.append({"label": row["label"], "arm": row["arm"],
                        "run_id": row["run_id"], "report_sha256": sha(report_raw),
                        "recovered_scalar": scalar, "sage_scalar_replay": True})
    assert len(checked) == 48
    receipt = {
        "kind": "s3_n41_c1200_c1800_confirmation24_independent_sage_scalar_replay_v1",
        "freeze_sha256": sha(freeze_path.read_bytes()),
        "ledger_sha256": sha(rows_raw),
        "replay_script_sha256": sha(Path(__file__).read_bytes()),
        "sage_runtime_info_sha256": sha(runtime_path.read_bytes()),
        "sage_version": SAGE_VERSION,
        "run_count": len(checked), "target_count": len(target_map),
        "all_scalars_independently_replayed": True,
        "runs": checked,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path),
                      "runs": len(checked), "targets": len(target_map)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger.resolve(), args.runtime_info.resolve(), args.receipt.resolve())
