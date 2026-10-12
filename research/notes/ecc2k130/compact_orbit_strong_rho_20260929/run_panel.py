#!/usr/bin/env python3
"""Four-arm cold point-only panel with balanced order and retained failures."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"))
from run_panel import child, jsonlines, sha  # noqa: E402


def run(n: int, length: int, rho_v2: Path, rho_v3: Path, ic: Path,
        out: Path) -> None:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    order = tuple(tuple(block) for block in frozen["arm_order"])
    assert frozen["evaluation_repetitions"] == len(order)
    assert all(set(block) == {"ic", "rho_v2", "rho_poly", "rho_normal"}
               for block in order)
    assert frozen["per_arm_timeout_seconds"] == 900
    assert frozen["per_arm_address_space_limit_bytes"] == 5 * 1024**3
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    points = HERE / spec["points_file"]
    fixture = HERE / spec["fixture_file"]
    assert sha(points) == spec["points_sha256"]
    assert sha(fixture) == spec["fixture_sha256"]
    assert len(jsonlines(fixture)) == length
    assert not out.exists(), "never overwrite a run directory"
    out.mkdir(parents=True)
    os.environ.pop("KIC_RHO_CANON_BACKEND", None)
    cpuinfo_path = Path("/proc/cpuinfo")
    cpuinfo = cpuinfo_path.read_text() if cpuinfo_path.exists() else ""
    cpu_fields = {}
    for line in cpuinfo.splitlines():
        if ":" in line:
            name, value = line.split(":", 1)
            cpu_fields.setdefault(name.strip(), value.strip())
    meminfo_path = Path("/proc/meminfo")
    meminfo = meminfo_path.read_text() if meminfo_path.exists() else ""
    host = {
        "platform": platform.platform(), "uname": tuple(platform.uname()),
        "python": sys.version, "cpu_count": os.cpu_count(),
        "cpu_model": cpu_fields.get("model name", platform.processor()),
        "cpu_flags": cpu_fields.get("flags", "").split(),
        "mem_total_kib": next((int(line.split()[1]) for line in meminfo.splitlines()
                               if line.startswith("MemTotal:")), None),
        "rustflags": os.environ.get("RUSTFLAGS"),
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "rho_v2_binary_sha256": sha(rho_v2),
        "rho_v3_binary_sha256": sha(rho_v3),
        "ic_binary_sha256": sha(ic),
        "load_at_start": os.getloadavg(), "pid": os.getpid(),
    }
    (out / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    k = frozen["compact_k"][str(n)][str(length)]
    (out / "chosen_k.json").write_text(json.dumps(
        {"n": n, "L": length, "k": k,
         "selection": "prior_disjoint_tune_no_eval_retuning"}, sort_keys=True) + "\n")
    runs = []
    for block, block_order in enumerate(order):
        for arm in block_order:
            prefix = out / f"b{block}_{arm}"
            if arm == "ic":
                env = {
                    "KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                    "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                }
                command = [str(ic), f"construct:{n}:0:{k}", str(points),
                           "7", str(prefix.with_suffix(".target.jsonl"))]
            else:
                env = {
                    "KIC_RHO_POINT_INPUT": str(points),
                    "KIC_RHO_BATCH_CORPUS": spec["corpus"],
                    "KIC_RHO_DP_BITS": str(frozen["rho_dp_bits"]),
                }
                if arm != "rho_v2":
                    env["KIC_RHO_CANON_BACKEND"] = (
                        "poly_xfirst" if arm == "rho_poly" else "normal_basis")
                binary = rho_v2 if arm == "rho_v2" else rho_v3
                command = [str(binary), str(n), "0", "signed_frobenius",
                           str(length), str(spec["seed"])]
            item = child(command, env, prefix)
            item.update({"block": block, "arm": arm})
            runs.append(item)
            (out / "runs.json").write_text(json.dumps(runs, indent=2, sort_keys=True) + "\n")
    (out / "host.json").write_text(json.dumps(
        dict(host, load_at_end=os.getloadavg()), indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(37, 41, 53), required=True)
    parser.add_argument("--L", type=int, choices=(1, 1024), required=True)
    parser.add_argument("--rho-v2", type=Path, required=True)
    parser.add_argument("--rho-v3", type=Path, required=True)
    parser.add_argument("--ic", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.n, args.L, args.rho_v2.resolve(), args.rho_v3.resolve(),
        args.ic.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
