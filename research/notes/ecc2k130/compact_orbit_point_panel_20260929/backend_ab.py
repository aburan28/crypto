#!/usr/bin/env python3
"""Same-runner paired rho backend comparison on frozen public Q files."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from run_panel import HERE, ROOT, child, jsonlines, sha


def run(n: int, length: int, baseline: Path, candidate: Path, out: Path) -> None:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert {"n": n, "L": length} in frozen["backend_ab_cases"]
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    points = HERE / spec["points_file"]
    assert sha(points) == spec["points_sha256"]
    assert not out.exists()
    out.mkdir(parents=True)
    cpuinfo = Path("/proc/cpuinfo").read_text() if Path("/proc/cpuinfo").exists() else ""
    model = next((line.split(":", 1)[1].strip() for line in cpuinfo.splitlines()
                  if line.startswith("model name")), platform.processor())
    flags = next((line.split(":", 1)[1].strip().split()
                  for line in cpuinfo.splitlines() if line.startswith("flags")), [])
    host = {
        "platform": platform.platform(), "uname": tuple(platform.uname()),
        "cpu_model": model, "cpu_flags": flags, "cpu_count": os.cpu_count(),
        "python": sys.version,
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "baseline_binary_sha256": sha(baseline),
        "candidate_binary_sha256": sha(candidate),
        "load_at_start": os.getloadavg(),
    }
    (out / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    runs = []
    for block in range(frozen["backend_ab_repetitions"]):
        for arm in (("baseline", "candidate") if block % 2 == 0
                    else ("candidate", "baseline")):
            prefix = out / f"b{block}_{arm}"
            env = {
                "KIC_RHO_POINT_INPUT": str(points),
                "KIC_RHO_BATCH_CORPUS": spec["corpus"],
                "KIC_RHO_DP_BITS": str(frozen["rho_dp_bits"]),
            }
            binary = baseline if arm == "baseline" else candidate
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
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.n, args.L, args.baseline.resolve(), args.candidate.resolve(),
        args.out.resolve())


if __name__ == "__main__":
    main()
