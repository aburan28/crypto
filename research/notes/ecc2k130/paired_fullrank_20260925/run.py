#!/usr/bin/env python3
"""Run one preregistered cold same-Q IC/rho pair; preserve failures and replay."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

MAIN_REF = "3bec69ca636cebc71f828a3035c49b8258e13f5f"
SEEDS = {
    37: (202609250037, 202609250137, 202609250237),
    41: (202609250041, 202609250141, 202609250241),
}
ENV = {
    "KIC_RANK_SURPLUS": "0",
    "KIC_INCREMENTAL_RANK_CROSSCHECK": "1",
    "RAYON_NUM_THREADS": "1",
}
SOURCE_PATHS = (
    "Cargo.toml",
    "examples/koblitz_rank_fixture.rs",
    "examples/koblitz_rho_fixture.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/binary_ecc/curve.rs",
    "src/binary_ecc/f2m.rs",
    "src/cryptanalysis/koblitz_fast.rs",
    "src/cryptanalysis/koblitz_groebner.rs",
    "src/cryptanalysis/inherited_f4.rs",
    "src/cryptanalysis/gf2_elim.rs",
    "src/cryptanalysis/fx_hash.rs",
    "src/cryptanalysis/mod.rs",
    "docs/ic/calibration.json",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_one(name: str, cmd: list[str], root: Path, out: Path) -> dict:
    stdout_path, stderr_path = out / f"{name}.jsonl", out / f"{name}.stderr.txt"
    process_env = os.environ.copy()
    process_env.update(ENV)
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        started = time.monotonic_ns()
        process = subprocess.Popen(cmd, cwd=root, env=process_env, stdout=stdout, stderr=stderr)
        timed_out = False
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                break
            if (time.monotonic_ns() - started) / 1e9 > 120:
                process.kill()
                _, status, usage = os.wait4(process.pid, 0)
                timed_out = True
                break
            time.sleep(0.01)
        ended = time.monotonic_ns()
        code = os.waitstatus_to_exitcode(status)
        process.returncode = code
    return {
        "arm": name,
        "argv": cmd,
        "env": ENV,
        "returncode": code,
        "timed_out": timed_out,
        "wall_ms": (ended - started) / 1e6,
        "cpu_s": usage.ru_utime + usage.ru_stime,
        "peak_rss_bytes": usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024,
        "stdout_sha256": sha(stdout_path),
        "stderr_sha256": sha(stderr_path),
        "binary_sha256": sha(Path(cmd[0])),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n", type=int, choices=SEEDS)
    parser.add_argument("--seed-index", type=int, choices=(0, 1, 2), required=True)
    args = parser.parse_args()
    assert args.n in SEEDS
    root = args.root.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    seed = SEEDS[args.n][args.seed_index]
    commands = {
        "ic": [
            str(root / "target/release/examples/koblitz_rank_fixture"),
            str(args.n), "0", "1", "2", "13737",
            "signed_expanded", "independent", "pair_pair_16", "1", f"hash:{seed}",
        ],
        "rho": [
            str(root / "target/release/examples/koblitz_rho_fixture"),
            str(args.n), "0", "signed_frobenius", "1", "packed", "13737", f"hash:{seed}",
        ],
    }
    order = ("ic", "rho") if args.seed_index % 2 == 0 else ("rho", "ic")
    report = {
        "schema_version": "1",
        "pinned_main_ref": MAIN_REF,
        "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "source_note": "Current-main producer and dependency overlay on a 74cf8424 checkout; rerun on an unmodified clean main checkout before making performance claims.",
        "source_sha256": {path: sha(root / path) for path in SOURCE_PATHS},
        "host": platform.platform(),
        "python": sys.version,
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "n": args.n,
        "a": 0,
        "public_hash_seed": seed,
        "seed_index": args.seed_index,
        "order": order,
        "arms": [],
    }
    for arm in order:
        item = run_one(arm, commands[arm], root, out)
        report["arms"].append(item)
        print(arm, item["returncode"], item["wall_ms"], item["cpu_s"], item["peak_rss_bytes"], flush=True)
    manifest_path = out / "manifest.json"
    manifest_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if all(item["returncode"] == 0 and not item["timed_out"] for item in report["arms"]):
        verifier = Path(__file__).with_name("verify.py")
        with (out / "replay.stderr.txt").open("wb") as stderr:
            completed = subprocess.run(
                [sys.executable, str(verifier), str(out / "ic.jsonl"), str(out / "rho.jsonl"),
                 "--out", str(out / "replay.json")],
                cwd=root, stderr=stderr, stdout=subprocess.PIPE, text=True,
            )
        report["replay_returncode"] = completed.returncode
        if completed.returncode == 0:
            report["replay_verdict"] = json.loads((out / "replay.json").read_text())["verdict"]
        else:
            report["replay_verdict"] = "FAIL"
    else:
        report["replay_verdict"] = "NOT_RUN_PRODUCER_FAILURE"
    manifest_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    assert report["replay_verdict"] == "PASS", report["replay_verdict"]


if __name__ == "__main__":
    main()
