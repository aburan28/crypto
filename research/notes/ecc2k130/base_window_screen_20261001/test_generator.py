#!/usr/bin/env python3
"""Frozen-source toy and archived n41 controls for the window generator."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile

from verify_generator import ARCHIVE_SHA256, FROZEN_COMPACT_SHA256, ROOT, sha, verify

HERE = Path(__file__).resolve().parent
TOY_BASE = HERE / "development/n13_k3_frozen_construct.base.jsonl"
TOY_BASE_SHA256 = "69fa0b23ab0604ecc7259dbe88770bf0514bed26168cb16eaf9b86d23fa41545"
FROZEN_FAST_ARITH_SHA256 = "2a5c54bdedb6a0b1ffcf95dd7badd32c41d1482b53055e31fb9f28d514ac8c9e"


def invoke(generator: Path, directory: Path, name: str, n: int, columns: int,
           cap: int) -> tuple[subprocess.CompletedProcess[str], Path, Path, Path]:
    header = directory / f"{name}.header.jsonl"
    scan = directory / f"{name}.scan.jsonl"
    receipt = directory / f"{name}.receipt.json"
    command = [
        str(generator), str(n), "0", str(columns), "0", str(cap),
        str(header), str(scan), str(receipt),
    ]
    return subprocess.run(command, capture_output=True, text=True, check=False), header, scan, receipt


def test(generator: Path) -> dict:
    assert generator.is_file()
    assert sha(TOY_BASE) == TOY_BASE_SHA256
    assert sha(ROOT / "src/cryptanalysis/koblitz_fast_arith.rs") == FROZEN_FAST_ARITH_SHA256
    with tempfile.TemporaryDirectory(prefix="base-window-control-") as tmp:
        directory = Path(tmp)
        toy, toy_header, toy_scan, toy_receipt = invoke(generator, directory, "toy", 13, 3, 1000)
        assert toy.returncode == 0, toy.stderr
        toy_check = verify(toy_header, toy_scan, toy_receipt)
        assert toy_header.read_bytes() == TOY_BASE.read_bytes()
        assert toy_check["orbit_members_checked"] == 78

        n41, n41_header, n41_scan, n41_receipt = invoke(generator, directory, "n41", 41, 255, 100000)
        assert n41.returncode == 0, n41.stderr
        n41_check = verify(n41_header, n41_scan, n41_receipt,
                           archived_control=True, protocol=True)
        assert n41_check["header_sha256"] == "33dd5e81eadcd773db8c3126c356ba6011c0fa7a02dfad62b411ce28e8a791c3"
        assert n41_check["archived_v2_compact_source_sha256"] == FROZEN_COMPACT_SHA256
        assert n41_check["accepted_orbits_replayed"] == 255
        assert n41_check["orbit_members_checked"] == 20910

        cap, cap_header, cap_scan, cap_receipt = invoke(generator, directory, "cap", 41, 255, 1)
        assert cap.returncode == 1 and not cap_header.exists()
        cap_check = verify(None, cap_scan, cap_receipt)
        assert cap_check["generator_status"] == "raw_x_cap"
        before = (sha(cap_scan), sha(cap_receipt))
        repeat, *_ = invoke(generator, directory, "cap", 41, 255, 1)
        assert repeat.returncode == 1 and "must not already exist" in repeat.stderr
        assert before == (sha(cap_scan), sha(cap_receipt))

        corrupted = directory / "corrupted.scan.jsonl"
        rows = [json.loads(line) for line in n41_scan.read_text().splitlines()]
        selected = next(row for row in rows if row.get("selected") is True)
        selected["selected"] = False
        corrupted.write_text("".join(json.dumps(row) + "\n" for row in rows))
        try:
            verify(n41_header, corrupted, n41_receipt, protocol=True)
        except AssertionError:
            mutation_rejected = True
        else:
            mutation_rejected = False
        assert mutation_rejected

        return {
            "schema": "koblitz-base-window-source-control-v1",
            "status": "PASS",
            "classification": "source_parity_only_no_scored_q_or_window_comparison",
            "frozen_compact_source_sha256": FROZEN_COMPACT_SHA256,
            "frozen_fast_arith_sha256": FROZEN_FAST_ARITH_SHA256,
            "toy": {
                "n": 13, "columns": 3,
                "frozen_construct_header_sha256": TOY_BASE_SHA256,
                "generated_header_sha256": toy_check["header_sha256"],
                "scan_decisions_replayed": toy_check["scan_decisions_replayed"],
                "orbit_members_checked": toy_check["orbit_members_checked"],
            },
            "n41": {
                "n": 41, "columns": 255,
                "archived_v2_tar_sha256": ARCHIVE_SHA256,
                "archived_v2_header_sha256": n41_check["archived_v2_header_sha256"],
                "generated_header_sha256": n41_check["header_sha256"],
                "scan_decisions_replayed": n41_check["scan_decisions_replayed"],
                "orbit_members_checked": n41_check["orbit_members_checked"],
            },
            "raw_x_cap_failure_replayed": cap_check["generator_status"] == "raw_x_cap",
            "existing_output_rejected_without_change": True,
            "scan_mutation_rejected": mutation_rejected,
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = test(args.generator.resolve())
    data = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.out:
        args.out.write_text(data)
    print(data, end="")


if __name__ == "__main__":
    main()
