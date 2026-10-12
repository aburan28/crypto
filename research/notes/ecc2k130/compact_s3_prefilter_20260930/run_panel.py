#!/usr/bin/env python3
"""Five cold, same-Q S3 blocked-filter/off/rho blocks with full receipts."""
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


ARMS = ("off_a", "filter", "off_b", "rho_normal")


def run(n: int, rho: Path, batch: Path, out: Path, cpu: int | None) -> None:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    source_lock_path = HERE / "SOURCE_FROZEN.json"
    assert sha(source_lock_path) == frozen["source_lock_sha256"]
    assert frozen["schema"] == "compact-s3-prefilter-evaluation-freeze-v1"
    for key, value in json.loads(source_lock_path.read_text()).items():
        if key != "schema":
            assert frozen[key] == value, key
    assert frozen["per_arm_timeout_seconds"] == 900
    assert frozen["per_arm_address_space_limit_bytes"] == 5 * 1024**3
    assert frozen["evaluation_repetitions"] == 5
    assert frozen["arm_order_offsets"] == [0, 1, 2, 3, 0]
    for filename, digest in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == digest, filename
    spec = frozen["specs"][f"n{n}_L1024_eval"]
    points = HERE / spec["points_file"]
    fixture = HERE / spec["fixture_file"]
    assert sha(points) == spec["points_sha256"]
    assert sha(fixture) == spec["fixture_sha256"]
    assert len(jsonlines(fixture)) == spec["L"] == 1024
    assert not out.exists(), "never overwrite an evaluation run"
    out.mkdir(parents=True)
    os.environ.pop("KIC_RHO_CANON_BACKEND", None)
    os.environ.pop("KIC_S3_BATCH_WINDOW", None)
    os.environ.pop("KIC_S3_PREFILTER", None)
    if sys.platform == "linux":
        assert cpu is not None and cpu in os.sched_getaffinity(0)
    cpuinfo = Path("/proc/cpuinfo")
    fields = {}
    for line in (cpuinfo.read_text() if cpuinfo.exists() else "").splitlines():
        if ":" in line:
            name, value = line.split(":", 1)
            fields.setdefault(name.strip(), value.strip())
    meminfo = Path("/proc/meminfo")
    memory = meminfo.read_text() if meminfo.exists() else ""
    host = {
        "platform": platform.platform(), "uname": tuple(platform.uname()),
        "python": sys.version, "cpu_count": os.cpu_count(),
        "cpu_model": fields.get("model name", platform.processor()),
        "cpu_flags": fields.get("flags", "").split(),
        "mem_total_kib": next((int(line.split()[1]) for line in memory.splitlines()
                               if line.startswith("MemTotal:")), None),
        "rustflags": os.environ.get("RUSTFLAGS"),
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "rho_binary_sha256": sha(rho),
        "batch_binary_sha256": sha(batch),
        "reserved_cpu": cpu,
        "load_at_start": os.getloadavg(), "pid": os.getpid(),
    }
    (out / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    k = frozen["k_by_n"][str(n)]
    (out / "grid.json").write_text(json.dumps({"n": n, "L": 1024,
        "k": k, "selection": "pre_registered_W64_incumbent_on_new_disjoint_corpus",
        "policies": list(ARMS)},
        indent=2, sort_keys=True) + "\n")
    runs = []
    for block, offset in enumerate(frozen["arm_order_offsets"]):
        order = ARMS[offset:] + ARMS[:offset]
        for arm in order:
            prefix = out / f"b{block}_{arm}"
            if arm != "rho_normal":
                assert arm in ("off_a", "filter", "off_b")
                env = {
                    "KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                    "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                    "KIC_S3_BATCH_WINDOW": "64",
                    "KIC_S3_PREFILTER": "blocked" if arm == "filter" else "off",
                }
                command = [str(batch), f"construct:{n}:0:{k}", str(points),
                           "7", str(prefix.with_suffix(".target.jsonl"))]
            else:
                assert arm == "rho_normal"
                env = {
                    "KIC_RHO_POINT_INPUT": str(points),
                    "KIC_RHO_BATCH_CORPUS": spec["corpus"],
                    "KIC_RHO_DP_BITS": str(frozen["rho_dp_bits"]),
                    "KIC_RHO_CANON_BACKEND": "normal_basis",
                }
                command = [str(rho), str(n), "0", "signed_frobenius",
                           "1024", str(spec["seed"])]
            if cpu is not None:
                command = ["taskset", "-c", str(cpu), *command]
            item = child(command, env, prefix)
            item.update({"block": block, "arm": arm})
            runs.append(item)
            (out / "runs.json").write_text(json.dumps(runs, indent=2, sort_keys=True) + "\n")
    (out / "host.json").write_text(json.dumps(
        dict(host, load_at_end=os.getloadavg()), indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(41, 53), required=True)
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cpu", type=int)
    args = parser.parse_args()
    run(args.n, args.rho.resolve(), args.batch.resolve(), args.out.resolve(), args.cpu)


if __name__ == "__main__":
    main()
