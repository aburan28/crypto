#!/usr/bin/env python3
"""Regression on the six immutable v1 pin-order counterexamples."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tarfile
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"))
from verify_panel import check_target as old_check_target  # noqa: E402
from check_pairing import check_target  # noqa: E402

RAW = ROOT / ("research/notes/ecc2k130/disjoint_cold_outcome_20261001/"
              "evidence_run_36794339148/raw/n37_L1024.tar.gz")
RAW_SHA = "f70d48bc110eec0514e1913c2657d8ed8221d4a83b70635c69da546659abca5d"
OLD_FIXTURES = ROOT / "research/notes/ecc2k130/disjoint_cold_q_20260930/fixtures"
CASES = {2: (541, 1006), 3: (435,)}


def member_rows(bundle: tarfile.TarFile, name: str) -> list[dict]:
    member = bundle.extractfile(f"disjoint-cold-n37_L1024/n37_L1024/{name}")
    assert member is not None
    return [json.loads(line) for line in member.read().splitlines() if line.strip()]


class PairingRegression(unittest.TestCase):
    def test_archived_mismatches_and_bad_pin(self) -> None:
        self.assertEqual(hashlib.sha256(RAW.read_bytes()).hexdigest(), RAW_SHA)
        checked = 0
        with tarfile.open(RAW, "r:gz") as bundle:
            for block, indices in CASES.items():
                fixture_path = OLD_FIXTURES / f"n37_L1024_b{block:02d}.fixture.jsonl"
                fixtures = [json.loads(line) for line in fixture_path.read_bytes().splitlines()]
                for arm in ("ic_a", "ic_b"):
                    prefix = f"b{block:02d}_{arm}"
                    base, = member_rows(bundle, f"{prefix}.base.jsonl")
                    rank = member_rows(bundle, f"{prefix}.rank.jsonl")
                    targets = member_rows(bundle, f"{prefix}.target.jsonl")
                    curve = Curve(Field(37, base["field_modulus_low_terms"]), 0)
                    generator = tuple(fixtures[0]["generator"])
                    logs = rank[-1]["logs"]
                    for index in indices:
                        target = targets[index]
                        fixture = fixtures[index]
                        with self.assertRaises(AssertionError):
                            old_check_target(target, fixture, base, logs, curve,
                                             generator, set())
                        check_target(target, fixture, base, logs, curve,
                                     generator, set())
                        wrong = dict(target)
                        wrong["pinned_intermediates"] = [
                            target["pinned_intermediates"][0] ^ 1,
                            target["pinned_intermediates"][1] ^ 1]
                        with self.assertRaises(AssertionError):
                            check_target(wrong, fixture, base, logs, curve,
                                         generator, set())
                        checked += 1
        self.assertEqual(checked, 6)


if __name__ == "__main__":
    unittest.main()
