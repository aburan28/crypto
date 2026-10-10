"""Independently replay the C1200 T220 direct and prior prefetch correctness diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "koblitz-s3-c1200-fold-n41-confirmation48-20261010"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main(direct_path: Path, prefetch_path: Path, runtime_path: Path,
         receipt_path: Path) -> None:
    assert not receipt_path.exists()
    runtime = json.loads(runtime_path.read_bytes())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    parent = json.loads((PARENT / "FREEZE.json").read_bytes())
    target = next(item for item in parent["selected_targets"] if item["label"] == "T220")
    ring = PolynomialRing(GF(2), "t")
    t = ring.gen()
    field = GF(2**41, name="z", modulus=t**41 + t**3 + 1)
    powers = [field.gen()**bit for bit in range(41)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << 41
        out = field(0)
        while word:
            bit = (word & -word).bit_length() - 1
            out += powers[bit]
            word &= word - 1
        return out

    curve = EllipticCurve(field, [1, 0, 0, 0, 1])

    def point(words: list[int]):
        assert len(words) == 2
        return curve(from_word(words[0]), from_word(words[1]))

    manifest = json.loads((PARENT / parent["arms"]["expanded"]["candidate_manifest_relative"]).read_bytes())
    generator = point(manifest["curve"]["generator_G"])
    r = manifest["curve"]["subgroup_order_r"]
    assert r == 549756390943 and generator != curve(0) and r * generator == curve(0)
    checked = {}
    for arm, path in (("direct", direct_path), ("prefetch", prefetch_path)):
        raw = path.read_bytes()
        report = json.loads(raw)
        assert report["n"] == 41 and report["a"] == 0
        assert report["orbit_columns"] == 1200 and report["prepared_rank"] == 1200
        assert report["factor_base_points"] == 98400
        assert report["factor_base_digest"] == parent["arms"]["expanded"]["factor_base_digest"]
        assert report["root_table_slots"] == parent["arms"]["expanded"]["root_table_slots"]
        assert report["group_verified"] is True
        assert report["target"] == target["public_point"]
        scalar = report["recovered_scalar"]
        assert scalar == target["fixture_scalar"] and 0 < scalar < r
        assert scalar * generator == point(target["public_point"])
        if arm == "direct":
            assert report["target_lookup_policy"] == "immediate_root_order_without_prefetch"
        checked[arm] = {
            "report_sha256": sha(raw), "recovered_scalar": scalar,
            "s3_calls": report["target_s3_calls"],
            "online_ms_diagnostic": report["timing_ms"]["target_online_after_reusable_setup"],
            "pdp_ms_diagnostic": report["timing_ms"]["target_pdp_charged"],
            "sage_scalar_replay": True,
        }
    receipt = {
        "kind": "s3_n41_c1200_directlookup_t220_sage_correctness_v1",
        "parent_freeze_sha256": sha((PARENT / "FREEZE.json").read_bytes()),
        "runtime_info_sha256": sha(runtime_path.read_bytes()),
        "replay_script_sha256": sha(Path(__file__).read_bytes()),
        "sage_version": SAGE_VERSION,
        "workload_id": target["workload_id"],
        "target": target["public_point"],
        "checked": checked,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"direct": checked["direct"]["sage_scalar_replay"],
                      "prefetch": checked["prefetch"]["sage_scalar_replay"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct", type=Path, required=True)
    parser.add_argument("--prefetch", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.direct.resolve(), args.prefetch.resolve(),
         args.runtime_info.resolve(), args.receipt.resolve())
