"""Independently replay every ARM PMULL fold screen scalar in Sage."""

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


def main(ledger: Path, runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_bytes())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    ledger_raw = ledger.read_bytes()
    rows = [json.loads(line) for line in ledger_raw.splitlines()]
    assert len(rows) == freeze["planned_run_count"] == 32
    assert [(r["n"], r["label"], r["arm"], r["repetition"]) for r in rows] == [
        (r["n"], r["label"], r["arm"], r["repetition"]) for r in freeze["run_order"]]
    targets = {(t["n"], t["label"]): t for t in freeze["selected_targets"]}
    ring = PolynomialRing(GF(2), "t")
    t = ring.gen()
    curves = {}
    for n, modulus in ((41, t**41 + t**3 + 1),
                       (53, t**53 + t**6 + t**2 + t + 1)):
        field = GF(2**n, name=f"z{n}", modulus=modulus)
        powers = [field.gen()**bit for bit in range(n)]

        def from_word(word: int, *, degree=n, f=field, p=powers):
            assert isinstance(word, int) and 0 <= word < 1 << degree
            value = f(0)
            while word:
                bit = (word & -word).bit_length() - 1
                value += p[bit]
                word &= word - 1
            return value

        curve = EllipticCurve(field, [1, 0, 0, 0, 1])

        def point(words: list[int], *, c=curve, decode=from_word):
            assert len(words) == 2
            return c(decode(words[0]), decode(words[1]))

        parent = (HERE.parent / "koblitz-s3-targetkernel-n41-20261010" if n == 41
                  else HERE.parent / "koblitz-s3-c1800-matchedsource-n53-20261010")
        old_arm = "improved" if n == 41 else "expanded"
        parent_freeze = json.loads((parent / "FREEZE.json").read_bytes())
        manifest_path = parent / parent_freeze["arms"][old_arm]["candidate_manifest_relative"]
        record = json.loads(manifest_path.read_bytes())
        generator = point(record["curve"]["generator_G"])
        order = record["curve"]["subgroup_order_r"]
        assert generator != curve(0) and order * generator == curve(0)
        curves[n] = (point, generator, order)
    checked = []
    for row in rows:
        n, label, arm = row["n"], row["label"], row["arm"]
        target = targets[(n, label)]
        frozen = freeze["curves"][str(n)]["arms"][arm]
        stem = f"N{n}_{label}_{arm}_R1"
        report_raw = (ledger.parent / f"{stem}.json").read_bytes()
        report = json.loads(report_raw)
        assert row["exit_code"] == 0 and row["verified_native_and_fixture"] is True
        assert row["report_sha256"] == sha(report_raw)
        assert row["candidate_id"] == frozen["candidate_id"]
        assert row["gf2_arm_pmull_fold"] == freeze["arms_env"][arm]
        assert row["workload_id"] == target["workload_id"]
        assert report["group_verified"] is True
        assert report["target"] == target["public_point"]
        assert report["factor_base_digest"] == freeze["curves"][str(n)]["factor_base_digest"]
        assert len(report["target_relation_indices"]) == 4
        point, generator, order = curves[n]
        scalar = int(report["recovered_scalar"])
        assert 0 < scalar < order
        assert scalar == target["fixture_scalar"]
        assert scalar * generator == point(target["public_point"])
        checked.append({"n": n, "label": label, "arm": arm,
                        "run_id": row["run_id"], "report_sha256": sha(report_raw),
                        "recovered_scalar": scalar, "sage_scalar_replay": True})
    receipt = {
        "kind": "s3_pmull_fold_screen_sage_replay_v1",
        "freeze_sha256": sha(freeze_path.read_bytes()),
        "ledger_sha256": sha(ledger_raw),
        "replay_script_sha256": sha(Path(__file__).read_bytes()),
        "sage_runtime_info_sha256": sha(runtime_path.read_bytes()),
        "sage_version": SAGE_VERSION,
        "run_count": len(checked), "target_count": len(targets),
        "all_scalars_independently_replayed": True,
        "runs": checked,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"runs": len(checked), "targets": len(targets)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger.resolve(), args.runtime_info.resolve(), args.receipt.resolve())
