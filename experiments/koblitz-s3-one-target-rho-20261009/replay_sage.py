"""Independently replay every paired public-point scalar in checked Sage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
CROSSCURVE = HERE.parent / "koblitz-s3-crosscurve-20261009"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(ledgers: list[Path], runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    freeze = json.loads((HERE / "FREEZE.json").read_text())
    plan_bytes = (CROSSCURVE / "FREEZE.json").read_bytes()
    assert sha(plan_bytes) == freeze["crosscurve_plan_sha256"]
    plan = json.loads(plan_bytes)
    targets = {(t["n"], t["label"]): t for t in plan["targets"]}
    ledger_bytes = [path.read_bytes() for path in ledgers]
    rows = [json.loads(line) for data in ledger_bytes for line in data.splitlines()]
    assert rows and len(rows) % 8 == 0
    checked = []
    for n in (41, 53):
        cfg = plan["curves"][str(n)]
        polynomial = PolynomialRing(GF(2), "t")
        t = polynomial.gen()
        field = GF(2**n, name="z", modulus=sum(t**p for p in cfg["modulus_powers"]))
        powers = [field.gen() ** bit for bit in range(n)]

        def from_word(word: int):
            assert 0 <= word < 1 << n
            result = field(0)
            while word:
                bit = (word & -word).bit_length() - 1
                result += powers[bit]
                word &= word - 1
            return result

        curve = EllipticCurve(field, [field(1), field(0), field(0), field(0), field(1)])
        generator = curve(*(from_word(word) for word in cfg["generator"]))
        assert cfg["r"] * generator == curve(0)
        for row in (item for item in rows if item["n"] == n):
            target = targets[(n, row["label"])]
            point = curve(*(from_word(word) for word in target["public_point"]))
            assert target["fixture_scalar"] * generator == point
            assert row["verified_native_and_fixture"] and row["phase_ok"]
            if row["arm"] == "s3":
                recovered = row["report"]["recovered_scalar"]
                assert row["report"]["target"] == target["public_point"]
            else:
                fixture = row["report"]["fixture"]
                recovered = fixture["recovered_fixture_scalar"]
                assert fixture["published_q"] == target["public_point"]
            assert recovered == target["fixture_scalar"]
            assert recovered * generator == point
            checked.append({"run_id": row["run_id"], "n": n,
                            "label": row["label"], "arm": row["arm"],
                            "independent_scalar_replay": True})
    assert len(checked) == len(rows)
    result = {"kind": "s3_rho_one_target_sage_replay_v1",
              "sage_version": SAGE_VERSION,
              "runtime_info_sha256": sha(runtime_path.read_bytes()),
              "crosscurve_plan_sha256": sha(plan_bytes),
              "ledger_sha256": [sha(data) for data in ledger_bytes],
              "run_count": len(rows), "target_count": len(rows) // 4,
              "all_scalars_independently_replayed": True,
              "runs": checked}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path), "runs": len(rows)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, nargs="+", required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.runtime_info, args.receipt)
