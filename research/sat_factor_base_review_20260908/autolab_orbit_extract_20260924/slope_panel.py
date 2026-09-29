#!/usr/bin/env python3
"""Compact-orbit extraction over the matched growing-n shared-log fixtures.

Runs `koblitz_s5_sat_instance` in `scalar` mode on every published rho fixture
scalar of the 2026-09-12 growing-n panel (n in {37,41,53}, 4 blocks x 32),
against the same retained factor base header. Public synthetic only.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PANEL = ROOT / "research/sat_factor_base_review_20260908/autolab_shared_log_scaling_20260912/runs"
OUT = Path(__file__).resolve().parent / "slope_panel"
BIN = ROOT / "target/release/examples/koblitz_s5_sat_instance"
ENV = {
    "KIC_ALGEBRA_ENCODING": "orbit_factorized",
    "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
    "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
    "KIC_ORBIT_REP_ENCODING": "one_hot",
}


def run_one(n, base_path, scalar, workdir):
    fixture_path = workdir / "fixture.json"
    env = dict(os.environ, **ENV, KIC_FACTOR_BASE_JSONL=str(base_path), KIC_FIXTURE_PATH=str(fixture_path))
    started = time.perf_counter()
    proc = subprocess.run(
        ["/usr/bin/time", "-l", str(BIN), str(n), "0", "1", "10", "scalar", str(scalar), "2000", "1", "internal"],
        capture_output=True, env=env,
    )
    wall_ms = (time.perf_counter() - started) * 1000.0
    match = re.search(rb"(\d+)\s+maximum resident set size", proc.stderr)
    peak_rss = int(match.group(1)) if match else None
    return proc.returncode, proc.stdout, proc.stderr, wall_ms, peak_rss, fixture_path


def main():
    if not BIN.exists():
        sys.exit(f"missing {BIN}; cargo build --release --example koblitz_s5_sat_instance")
    OUT.mkdir(parents=True, exist_ok=True)
    fields = [int(v) for v in (sys.argv[1:] or ["37", "41", "53"])]
    for n in fields:
        for block_dir in sorted((PANEL / f"n{n}").glob("block_*")):
            block = int(block_dir.name.split("_")[1])
            record_path = OUT / f"n{n}_block{block:02d}.jsonl"
            if record_path.exists():
                continue
            header = (block_dir / "direct/stdout.txt").open("rb").readline()
            base_hash = json.loads(header)["base_hash"]
            with tempfile.TemporaryDirectory() as tmp:
                workdir = Path(tmp)
                base_path = workdir / "base.jsonl"
                base_path.write_bytes(header)
                records = []
                for line in (block_dir / "rho/stdout.txt").read_text().splitlines():
                    rho = json.loads(line)
                    code, stdout, stderr, wall_ms, peak_rss, fixture_path = run_one(
                        n, base_path, rho["published_fixture_scalar"], workdir
                    )
                    rec = {
                        "n": n, "block": block, "fixture_index": rho["fixture_index"],
                        "published_fixture_scalar": rho["published_fixture_scalar"],
                        "published_q": rho["published_q"], "base_hash": base_hash,
                        "rho_total_ms": rho["total_ms"], "rho_walk_ms": rho["walk_ms"],
                        "exit_code": code, "process_wall_ms": wall_ms, "peak_rss_bytes": peak_rss,
                    }
                    if code == 0:
                        d = json.loads(stdout)
                        fx = json.loads(fixture_path.read_text())
                        e = d["compact_orbit_extraction"]
                        rec.update({
                            "generator": fx["generator"], "target": fx["target"],
                            "field_modulus_low_terms": fx["field_modulus_low_terms"],
                            "subgroup_order": fx["subgroup_order"],
                            "target_matches_published_q": fx["target"] == rho["published_q"],
                            "factor_base_input_hash": d["factor_base_input_hash"],
                            "sat_verdict": d["decomposition_verdict"], "sat_conflicts": d["conflicts"],
                            "valid_x_tuples": d["valid_x_tuples"], "invalid_group_lifts": d["invalid_group_lifts"],
                            "pair_table_entries": d["pair_table_entries"], "edge_selectors": e["edge_selectors"],
                            "group_valid": e["group_valid"], "extract_ms": e["extract_ms"],
                            "extract_trials": e["trials"], "index_entries": e["index_entries"],
                            "regular_states": d["lazy_relative_support"]["regular_states"],
                            "scan_ms": d["lazy_relative_support"]["scan_ms"],
                            "base_ms": d["base_ms"], "encoding_ms": d["encoding_ms"],
                            "solve_ms": d["solve_ms"], "x_codes": e["x_codes"],
                            "pinned_intermediates": e["pinned_intermediates"],
                        })
                    records.append(rec)
                    print(n, block, rec["fixture_index"], rec.get("sat_verdict"), round(rec.get("extract_ms", -1), 1), flush=True)
            record_path.write_text("".join(json.dumps(r) + "\n" for r in records))


if __name__ == "__main__":
    main()
