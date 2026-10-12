#!/usr/bin/env python3
"""One phase of one sweep cell (RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md).

usage: sweep_cell.py <phase> <n> <L> <K> <corpus>
  phase native         rho R3 (32 lanes) then IC, uninstrumented; writes
                       scalars.txt (the rho corpus' planted scalars, fed to IC)
  phase callgrind_rho  rho R3 under valgrind --tool=callgrind
  phase callgrind_ic   IC under valgrind --tool=callgrind

Everything lands in cell_n<n>_L<L>_K<K>/ next to this script.
"""
import json
import os
import resource
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
RHO = os.path.join(ROOT, "target", "release", "examples", "koblitz_rho_batch_ks_strong")
IC = os.path.join(ROOT, "target", "release", "examples", "koblitz_orbit_dlp_fast")
BATCH_SEED = "531310"


def uptime():
    return subprocess.run(["uptime"], capture_output=True, text=True).stdout.strip()


def timed(cmd, stdout_path, stderr_path, env=None):
    """Run cmd, return (rc, wall_s, child_maxrss_bytes)."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    t0 = time.perf_counter()
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        rc = subprocess.run(cmd, stdout=out, stderr=err, env=env).returncode
    wall = time.perf_counter() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # ru_maxrss is KiB on Linux and is the max over all children so far.
    return rc, wall, max(after, before) * 1024


def rho_env(corpus):
    env = dict(os.environ)
    env.update(
        KIC_RHO_RUNG="3",
        KIC_RHO_LANES="32",
        KIC_RHO_DP_BITS="4",
        KIC_RHO_BATCH_CORPUS=corpus,
    )
    return env


def rho_cmd(n, L):
    return [RHO, str(n), "0", "signed_frobenius", str(L), BATCH_SEED]


def ic_cmd(n, K, out_jsonl):
    return [IC, f"construct:{n}:0:{K}", "scalars.txt", "7", out_jsonl]


def parse_rho(path):
    fixtures, summary = [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            if obj.get("kind") == "rho_ks_batch_fixture":
                fixtures.append(obj)
            elif obj.get("kind") == "rho_ks_batch_summary":
                summary = obj
    return fixtures, summary


def callgrind_ir(stderr_path):
    with open(stderr_path) as fh:
        for line in fh:
            if "I   refs:" in line:
                return int(line.split("refs:")[1].strip().replace(",", ""))
    return None


def main():
    phase, n, L, K, corpus = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    cell = os.path.join(HERE, f"cell_n{n}_L{L}_K{K}")
    os.makedirs(cell, exist_ok=True)
    os.chdir(cell)
    result = {"phase": phase, "n": n, "L": L, "K": K, "corpus": corpus, "batch_seed": BATCH_SEED}
    result["load_before"] = uptime()

    if phase == "native":
        rc, wall, rss = timed(rho_cmd(n, L), "rho_native.jsonl", "rho_native.stderr.log", rho_env(corpus))
        fixtures, summary = parse_rho("rho_native.jsonl")
        with open("scalars.txt", "w") as fh:
            for f in fixtures:
                fh.write(f"{f['published_fixture_scalar']}\n")
        rho_ok = (
            rc == 0
            and len(fixtures) == L
            and summary is not None
            and summary["all_verified"] is True
            and all(f["verified"] and f["recovered_fixture_scalar"] == f["published_fixture_scalar"] for f in fixtures)
        )
        result["rho"] = {
            "rc": rc,
            "wall_s": wall,
            "child_maxrss_bytes": rss,
            "gates_G1": rho_ok,
            "total_walk_steps": summary and summary["total_walk_steps"],
            "table_entries": summary and summary["table_entries"],
            "in_process_ms": summary and summary["in_process_ms"],
        }
        rc, wall, rss = timed(ic_cmd(n, K, "ic_native.jsonl"), "ic_native.summary.json", "ic_native.stderr.log")
        try:
            ic = json.load(open("ic_native.summary.json"))
        except Exception:  # noqa: BLE001 - reported below
            ic = {}
        result["ic"] = {
            "rc": rc,
            "wall_s": wall,
            "gates_G1": rc == 0 and ic.get("targets_solved") == L and ic.get("targets_failed") == 0,
            "targets_solved": ic.get("targets_solved"),
            "targets_failed": ic.get("targets_failed"),
            "rank": ic.get("rank"),
            "orbit_columns": ic.get("orbit_columns"),
            "rank_failures": ic.get("rank_failures"),
            "peak_rss_bytes": ic.get("peak_rss_bytes"),
            "timing_ms": ic.get("timing_ms"),
            "root_table_entries": ic.get("root_table_entries"),
        }
    elif phase == "callgrind_rho":
        cmd = ["valgrind", "--tool=callgrind", "--callgrind-out-file=callgrind_rho.out"] + rho_cmd(n, L)
        rc, wall, rss = timed(cmd, "callgrind_rho.stdout.jsonl", "callgrind_rho.stderr.log", rho_env(corpus))
        fixtures, summary = parse_rho("callgrind_rho.stdout.jsonl")
        result["rho"] = {
            "rc": rc,
            "wall_s": wall,
            "ir": callgrind_ir("callgrind_rho.stderr.log"),
            "gates_G1": rc == 0 and len(fixtures) == L and summary is not None and summary["all_verified"] is True,
            "total_walk_steps": summary and summary["total_walk_steps"],
        }
    elif phase == "callgrind_ic":
        cmd = ["valgrind", "--tool=callgrind", "--callgrind-out-file=callgrind_ic.out"] + ic_cmd(n, K, "callgrind_ic.stdout.jsonl")
        rc, wall, rss = timed(cmd, "callgrind_ic.summary.json", "callgrind_ic.stderr.log")
        try:
            ic = json.load(open("callgrind_ic.summary.json"))
        except Exception:  # noqa: BLE001
            ic = {}
        result["ic"] = {
            "rc": rc,
            "wall_s": wall,
            "ir": callgrind_ir("callgrind_ic.stderr.log"),
            "gates_G1": rc == 0 and ic.get("targets_solved") == L and ic.get("targets_failed") == 0,
            "targets_solved": ic.get("targets_solved"),
            "targets_failed": ic.get("targets_failed"),
            "peak_rss_bytes": ic.get("peak_rss_bytes"),
        }
    else:
        raise SystemExit(f"unknown phase {phase}")

    result["load_after"] = uptime()
    with open(f"{phase}_result.json", "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k in ("phase", "n", "L", "K", "rho", "ic")}))


if __name__ == "__main__":
    main()
