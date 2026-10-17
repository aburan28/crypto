#!/usr/bin/env python3
"""Derive root22 configurations from independently replayed equal bases."""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINEGRID = ROOT / "research/n83_e2e_finegrid_20261016"
KNOWN_IDS = {
    1104: "IC1N83Ckb1a1fb183264PDP4rootRCguidedLAgaussTDdirectISO0hce7cffedf9c0",
    1120: "IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0h82b74c23b7bc",
    1107: "IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0hcb73613ebc5a",
    1112: "IC1N83Ckb1a1fb184592PDP4rootRCguidedLAgaussTDdirectISO0hbc5ec3072dc4",
}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finegrid = module("n83_root22_candidate_finegrid", FINEGRID / "materialize.py")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("finegrid_results", type=Path)
    args = parser.parse_args()
    results = args.finegrid_results.resolve()
    rows = json.loads((FINEGRID / "ic1_rows_38978215742.json").read_text())
    if (rows["selection_status"] != "SCREEN_COMPLETE"
            or rows["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or rows["raw_artifact_sha256"] !=
            "b07cef5779312fed6fb49c2af767d5d2c8cbdc1eca6c14de85ddd6753fb6290f"):
        raise ValueError("finegrid replay and provider binding differ")
    field, curve = finegrid.panel.curve_record()
    for k, expected_id in KNOWN_IDS.items():
        name = f"k{k}_a" if k in (1104, 1120) else f"k{k}"
        folder = results / "runs" / name
        summary = finegrid.panel.read(folder / "stdout.jsonl")
        dump = folder / "base_dump.jsonl"
        root11_id = rows["by_k"][str(k)]["candidate_id"]
        root11 = finegrid.panel.read(FINEGRID / f"candidate_{root11_id}.json")
        if (finegrid.sha256(dump) !=
                root11["factor_base"]["ordered_base_dump_sha256"]
                or summary["base_hash"] !=
                root11["factor_base"]["source_base_hash_blake3"]):
            raise ValueError(f"K{k}: finegrid base dump differs")
        summary["root_state_layout"] = "root_only83_22byte_states"
        summary["root_state_slot_bytes"] = 22
        summary["root_only11_enabled"] = False
        manifest = finegrid.dynamic.candidate(
            k, "root22_dynamic", summary, dump, field, curve)
        if manifest["candidate_id"] != expected_id:
            raise ValueError(f"K{k}: root22 identity differs")
        if k in (1107, 1112):
            finegrid.dynamic.write_once(HERE / f"candidate_{expected_id}.json",
                                        manifest)
        print(expected_id)


if __name__ == "__main__":
    main()
